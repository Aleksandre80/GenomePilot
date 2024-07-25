import pysam
import sys

input_bam_path = sys.argv[1]
output_bam_path = sys.argv[2]
chromosome = sys.argv[3]
position1 = int(sys.argv[4])  # Première position génomique, base 1
position2 = int(sys.argv[5])  # Seconde position génomique, base 1
min_quality = int(sys.argv[6])  # Qualité Phred minimale pour conserver un read
log_file_path = sys.argv[7]

input_bam = pysam.AlignmentFile(input_bam_path, "rb")
output_bam = pysam.AlignmentFile(output_bam_path, "wb", template=input_bam)

with open(log_file_path, 'w') as log_file:
    log_file.write("ReadName\tBaseAtPosition1\tQuality1\tBaseAtPosition2\tQuality2\tStatus\n")

    # Parcourir tous les reads qui couvrent au moins une des positions
    start = min(position1, position2) - 1
    end = max(position1, position2)
    for read in input_bam.fetch(chromosome, start, end):
        if not read.is_unmapped:
            qualities = []
            bases = ['N', 'N']
            positions = [position1, position2]
            status = "Kept"

            for i, pos in enumerate(positions):
                read_pos = None
                for qpos, rpos in read.get_aligned_pairs():
                    if rpos == pos - 1:
                        read_pos = qpos
                        break

                if read_pos is not None:
                    if read.query_qualities is not None and 0 <= read_pos < len(read.query_qualities):
                        qualities.append(read.query_qualities[read_pos])
                    else:
                        qualities.append(None)

                    if read.query_sequence is not None and 0 <= read_pos < len(read.query_sequence):
                        bases[i] = read.query_sequence[read_pos]
                else:
                    qualities.append(None)

            if any(q is not None and q < min_quality for q in qualities):
                status = "Removed"
            else:
                output_bam.write(read)

            log_file.write(f"{read.query_name}\t{bases[0]}\t{qualities[0]}\t{bases[1]}\t{qualities[1]}\t{status}\n")

input_bam.close()
output_bam.close()
