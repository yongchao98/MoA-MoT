import unicodedata

def reconstruct_sanskrit_ancestor():
  """
  This function provides the reconstructed Sanskrit ancestor in IPA.
  The reconstruction is based on the comparative method applied to the provided
  Hindi and Punjabi words.

  The logic is as follows:
  1. Hindi 'kãːs' and Punjabi 'kə̃d͡ʒʱ' are the descendants.
  2. The correspondence k:k -> *k.
  3. The correspondence H. ãːC : P. ə̃CC suggests a proto-form with *a + Nasal + Consonant Cluster.
  4. The correspondence H. s : P. d͡ʒʱ is key. This irregular correspondence points to a single
     Sanskrit cluster with divergent reflexes.
  5. The Punjabi 'd͡ʒʱ' regularly derives from Sanskrit '*dhy' (via MIA 'jjh').
  6. Let's posit the ancestor *kandhya.
  7. Path to Punjabi: Skt. *kandhya > MIA *kaṁjjha > Pbi. kə̃d͡ʒʱ. This is a regular development.
  8. Path to Hindi: Skt. *kandhya > MIA *kaṁjjha > standard H. *kā̃jh. The given form 'kãːs'
     reflects a known but irregular secondary change of jha > sa.
  9. This makes *kandhya the most probable reconstruction.
  """
  # The Sanskrit form is *kandhya.
  # In IPA, this is represented as:
  # k = k (voiceless velar stop)
  # a = a (open central unrounded vowel)
  # n = n̪ (dental nasal)
  # dh = d̪ʱ (voiced dental aspirated stop)
  # y = j (palatal approximant)
  # a = a (final vowel)
  ipa_reconstruction = "kand̪ʱja"
  
  # The final output requires printing the parts of the "equation"
  # Let's consider the equation to be the mapping from letters to IPA symbols.
  print("Reconstructed Sanskrit Form: *kandhya")
  print("IPA Breakdown of the reconstructed form:")
  print("k -> k")
  print("a -> a")
  print("n -> n̪")
  print("dh -> d̪ʱ")
  print("y -> j")
  print("a -> a")
  print("\nFinal Reconstructed Form in IPA:")
  print(ipa_reconstruction)

reconstruct_sanskrit_ancestor()