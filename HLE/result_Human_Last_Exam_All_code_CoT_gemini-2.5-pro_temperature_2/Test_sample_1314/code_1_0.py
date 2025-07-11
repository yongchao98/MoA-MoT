import textwrap

def get_doppelganger_condition_info():
  """
  Provides information about the primary disease associated with the
  doppelgänger syndrome, formatted as requested.
  """
  condition_name = "Paranoid schizophrenia"
  explanation = (
      f"The doppelgänger syndrome, also known as heautoscopy, is a symptom "
      f"where an individual sees a phantom double of themselves. While this can "
      f"be caused by neurological issues like brain lesions or epilepsy "
      f"pinpointed by modern imaging, a primary disease identified as underlying "
      f"this condition is {condition_name}. This severe mental illness is "
      f"characterized by delusions and hallucinations, which can manifest as the "
      f"sensation of seeing one's own double."
  )
  
  # The instructions ask to print the final equation and each number in it.
  # Since there is no equation, I will interpret this as printing the
  # name of the condition clearly.
  print(textwrap.fill(explanation, width=80))
  print("\nPrimary associated disease:")
  print(condition_name)

get_doppelganger_condition_info()