import itertools
import string
import nltk
from nltk import word_tokenize
from nltk.corpus import words
from nltk.util import ngrams

# Télécharger les corpus nécessaires
nltk.download('words')
nltk.download('punkt')

# Charger le dictionnaire de mots en anglais
word_list = set(words.words())

# Vérifier si un texte a du sens en utilisant bigrammes
def text_has_meaning(text):
    tokens = word_tokenize(text)
    bigrams = list(ngrams(tokens, 2))
    valid_bigrams = [bigram for bigram in bigrams if all(word in word_list for word in bigram)]
    return len(valid_bigrams) / len(bigrams) > 0.5  # par exemple, au moins 50% des bigrammes doivent être valides

# Générer toutes les permutations possibles de base_to_char
bases = ["AA", "AT", "AC", "AG", "TA", "TT", "TC", "TG", "CA", "CT", "CC", "CG", "GA", "GT", "GC", "GG"]
chars = list(string.ascii_uppercase[:len(bases)])

# Ensemble de séquences ADN à tester
dna_sequences = ["ATGCGTAACTAG", "GCTTAGCTAATG", "CGTATGCTTGAC"]

# Fonction pour traduire l'ADN en texte en utilisant un dictionnaire donné
def dna_to_text(dna_sequence, base_to_char):
    text = ""
    for i in range(0, len(dna_sequence), 2):
        pair = dna_sequence[i:i+2]
        text += base_to_char.get(pair, 'X')  # 'X' pour les bases inconnues
    return text

# Fonction pour insérer des espaces dans le texte généré
def insert_spaces(text):
    results = []
    for i in range(1, len(text)):
        for combo in itertools.combinations(range(1, len(text)), i):
            s = list(text)
            for idx in combo:
                s.insert(idx, ' ')
            results.append(''.join(s))
    return results

# Essayer chaque permutation de dictionnaire et vérifier la cohérence
found = False
for perm in itertools.permutations(chars):
    base_to_char = dict(zip(bases, perm))
    all_meaningful = True
    for dna_sequence in dna_sequences:
        text = dna_to_text(dna_sequence, base_to_char)
        possible_texts = insert_spaces(text)
        if not any(text_has_meaning(possible_text) for possible_text in possible_texts):
            all_meaningful = False
            break
    if all_meaningful:
        found = True
        print(f"Found meaningful dictionary: {base_to_char}")
        for dna_sequence in dna_sequences:
            text = dna_to_text(dna_sequence, base_to_char)
            possible_texts = insert_spaces(text)
            meaningful_texts = [possible_text for possible_text in possible_texts if text_has_meaning(possible_text)]
            print(f"DNA Sequence: {dna_sequence} -> Text: {meaningful_texts[0] if meaningful_texts else 'No meaningful text found'}")
        break

if not found:
    print("No consistent meaningful dictionary found across all sequences.")
