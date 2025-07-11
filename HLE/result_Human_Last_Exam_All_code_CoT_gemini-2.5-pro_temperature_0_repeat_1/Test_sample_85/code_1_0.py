import unicodedata

def get_reconstruction():
  """
  This function returns the reconstructed Sanskrit ancestor in IPA.
  The reconstruction is based on the analysis of the descendants:
  - Hindi: kãːs
  - Punjabi: kə̃d͡ʒʱ

  The reconstructed form is Skt. kaṅkṣa, which has different reflexes
  in the dialects leading to Hindi and Punjabi.
  - Path to Hindi: kaṅkṣa > MIA *kaṃsa > kãːs
  - Path to Punjabi: kaṅkṣa > MIA *kaṃjjha > kə̃d͡ʒʱ

  The IPA for kaṅkṣa is 'kaŋkʂa'.
  """
  # Sanskrit: k a ṅ k ṣ a
  # IPA:      k a ŋ k ʂ a
  reconstructed_ipa = "kaŋkʂa"
  return reconstructed_ipa

# Print the final reconstructed word in IPA
print(get_reconstruction())