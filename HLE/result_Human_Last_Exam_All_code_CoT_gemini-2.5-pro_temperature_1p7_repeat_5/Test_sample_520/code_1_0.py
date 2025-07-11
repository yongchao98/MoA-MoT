def find_profession():
  """
  This function provides information about Andre Vernault's true profession.
  """
  # Data about Andre Vernault based on historical records.
  # He was a Belgian refugee who arrived in the UK during WWII.
  # British intelligence suspected him of being a spy because he was evasive
  # about his background, fearing his real job would lead to ridicule.
  vernault_file = {
      "name": "Andre Vernault",
      "status": "Belgian refugee",
      "suspected_of": "Being a German spy",
      "concealed_true_profession": "Female impersonator"
  }

  # Retrieve the true profession from the data
  true_profession = vernault_file["concealed_true_profession"]

  # Print the answer
  print("During World War II, Andre Vernault concealed his true profession from interrogators.")
  print(f"His true profession was: {true_profession}")

find_profession()