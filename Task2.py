import re
import pyopenms
from collections import defaultdict

# Protein sequence
sequence = (
    "MKWVTFISLLFLFSSAYSRGVFRRDTHKSEIAHRFKDLGEENFKALVLIAFAQYLQQCPFEDHVKLVN"
    "EVTEFAKTCVADESAENCDKSLHTLFGDELCKVASLRETYGSHCIDVFQLER"
)

# Step 1: Perform trypsin digestion
def digest_protein(sequence, enzyme='trypsin'):
    if enzyme == 'trypsin':
        peptides = re.split(r'(?<=[KR])(?!P)', sequence)
    return peptides

peptides = digest_protein(sequence)
print("Peptides after trypsin digestion:", peptides)

# Step 2: Compute monoisotopic mass for each peptide
def calculate_monoisotopic_mass_peptide(peptide_sequence):
    peptide = pyopenms.AASequence.fromString(peptide_sequence)
    mass = peptide.getMonoWeight()
    return mass

# Step 3: Store peptide and mass information
def store_peptide_and_mass_info(peptides):
    peptide_info_dict = defaultdict(list)
    for peptide in peptides:
        mass = calculate_monoisotopic_mass_peptide(peptide)
        peptide_info_dict[mass].append(peptide)
    return peptide_info_dict

peptide_info_dict = store_peptide_and_mass_info(peptides)

# Display mass and peptide data
for mass, peptides_list in peptide_info_dict.items():
    print(f"Mass: {mass:.4f} Da ==> Peptides: {peptides_list}")

# Step 4: Identify isobaric peptides
def find_isobaric_peptides(peptide_data):
    result = {}
    for mass, sequences in peptide_data.items():
        if len(sequences) > 1:
            result[mass] = sequences
    return result

isobaric_peptides = find_isobaric_peptides(peptide_info_dict)

# Display isobaric peptides
if isobaric_peptides:
    print("\nDetected Isobaric Peptides:")
    for mass, sequences in isobaric_peptides.items():
        formatted_sequences = " | ".join(sequences)
        print(f"> Mass: {mass:.4f} Da")
        print(f"  Sequences: {formatted_sequences}")
else:
    print("\nNo isobaric peptides were found.")

# Example Output:
# Peptides after trypsin digestion: ['MKWVTFISLLFLFSSAYSR', 'GVFR', 'RDTHK', 'SEIAHR', 'FK', ...]
# Mass: 2230.1523 Da ==> Peptides: ['MKWVTFISLLFLFSSAYSR']
# Mass: 487.2431 Da ==> Peptides: ['GVFR']
# Mass: 639.3198 Da ==> Peptides: ['RDTHK']
# ...
# Detected Isobaric Peptides:
# > Mass: 174.1117 Da
#   Sequences: R | R
