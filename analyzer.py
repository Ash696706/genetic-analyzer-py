from Bio.Seq import Seq


# -------------------------------
# DNA ANALYSIS
# -------------------------------

def analyze_dna(dna):
    dna = Seq(dna.upper())

    length = len(dna)

    a_count = dna.count("A")
    t_count = dna.count("T")
    g_count = dna.count("G")
    c_count = dna.count("C")

    gc_content = ((g_count + c_count) / length) * 100

    return {
        "length": length,
        "A": a_count,
        "T": t_count,
        "G": g_count,
        "C": c_count,
        "GC%": round(gc_content, 2)
    }


# -------------------------------
# TRANSCRIPTION & TRANSLATION
# -------------------------------

def analyze_protein(dna):
    dna = Seq(dna.upper())

    rna = dna.transcribe()
    protein = rna.translate()

    stop_codons = protein.count("*")

    return rna, protein, stop_codons


# -------------------------------
# MUTATION ANALYSIS
# -------------------------------

def find_mutations(reference, sample):

    mutations = []

    if len(reference) != len(sample):
        print("Error: sequences must have the same length.")
        return mutations

    for position, (ref_base, sample_base) in enumerate(
        zip(reference, sample), start=1
    ):

        if ref_base != sample_base:

            mutations.append({
                "position": position,
                "reference": ref_base,
                "sample": sample_base
            })

    return mutations


# -------------------------------
# MAIN PROGRAM
# -------------------------------

reference = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"

sample = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAC"


# Analyze reference
reference_analysis = analyze_dna(reference)

# Analyze sample
sample_analysis = analyze_dna(sample)

# Protein analysis
reference_rna, reference_protein, reference_stops = analyze_protein(reference)

sample_rna, sample_protein, sample_stops = analyze_protein(sample)

# Mutation analysis
mutations = find_mutations(reference, sample)


# -------------------------------
# REPORT
# -------------------------------

print("=" * 50)
print("        🧬 GENETIC SEQUENCE ANALYZER")
print("=" * 50)

print("\nREFERENCE DNA")
print(reference)

print("\nSAMPLE DNA")
print(sample)


print("\n" + "-" * 50)
print("DNA ANALYSIS")
print("-" * 50)

print("DNA length:", reference_analysis["length"])

print(
    "Base counts:",
    "A =", reference_analysis["A"],
    "| T =", reference_analysis["T"],
    "| G =", reference_analysis["G"],
    "| C =", reference_analysis["C"]
)

print("GC content:", reference_analysis["GC%"], "%")


print("\n" + "-" * 50)
print("TRANSCRIPTION")
print("-" * 50)

print("mRNA:")
print(reference_rna)


print("\n" + "-" * 50)
print("TRANSLATION")
print("-" * 50)

print("Reference protein:")
print(reference_protein)

print("Sample protein:")
print(sample_protein)

print("Reference stop codons:", reference_stops)
print("Sample stop codons:", sample_stops)


print("\n" + "-" * 50)
print("MUTATION ANALYSIS")
print("-" * 50)

print("Total mutations:", len(mutations))

if len(mutations) == 0:

    print("No mutations detected.")

else:

    for mutation in mutations:

        print(
            f"Position {mutation['position']}: "
            f"{mutation['reference']} → {mutation['sample']}"
        )


print("\n" + "=" * 50)
print("Analysis complete.")
print("=" * 50)
