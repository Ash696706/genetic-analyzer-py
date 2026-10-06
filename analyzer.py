from Bio.Seq import Seq

# dna analysis
def analyze_dna(dna):
  dna = Seq(dna.upper())

  length = len(dna)

  a_count = dna.count("A")
  t_count = dna.count("T")
  g_count = dna.count("G")
  c_count = dna.count("C")

  gc_content = ((g_count + c_count)/length)*100

  return (
      "Length:", length,
      "A count:", a_count,
      "T count:", t_count,
      "G count:", g_count,
      "C count:", c_count,
      "GC% :", round(gc_content,2), "%"
  )
#dna sequence
print(analyze_dna("ATGCct"))


#transcription and translation
def analyze_dna(dna):
  dna = Seq(dna.upper())

  rna = dna.transcribe()
  protein = rna.translate()

  stop_codons = protein.count("*")

  return rna, protein, stop_codons
#dna sequence
print(analyze_dna("ATGCct"))
