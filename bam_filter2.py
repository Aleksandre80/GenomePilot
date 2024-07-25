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

    for read in input_bam.fetch():
        if not read.is_unmapped:
            positions = [position1, position2]
            pos_covered = [False, False]
            qualities = [None, None]
            bases = ['N', 'N']

            # Check coverage and quality at specified positions
            for i, pos in enumerate(positions):
                for qpos, rpos in read.get_aligned_pairs():
                    if rpos == pos - 1:
                        pos_covered[i] = True
                        if read.query_qualities is not None and 0 <= qpos < len(read.query_qualities):
                            qualities[i] = read.query_qualities[qpos]
                        if read.query_sequence is not None and 0 <= qpos < len(read.query_sequence):
                            bases[i] = read.query_sequence[qpos]
                        break

            # Determine read status based on quality checks
            if any((q is not None and q < min_quality) for q, covered in zip(qualities, pos_covered) if covered):
                status = "Removed"
            else:
                status = "Kept"
                output_bam.write(read)

            # Logging each read's handling
            log_file.write(f"{read.query_name}\t{bases[0]}\t{qualities[0]}\t{bases[1]}\t{qualities[1]}\t{status}\n")

input_bam.close()
output_bam.close()
