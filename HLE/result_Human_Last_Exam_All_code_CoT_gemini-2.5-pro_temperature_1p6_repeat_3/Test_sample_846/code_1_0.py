def explain_azobenzene_table():
  """
  This function describes the necessary behavior of a picnic table
  shaped like an azobenzene molecule to be functionally accurate.
  """
  sunrise_action = "At sunrise, when exposed to UV light from the sun, the table should change shape from its stable, extended 'trans' form to a bent 'cis' form."
  sunset_action = "At sunset, in the absence of UV light, the table should relax back from the bent 'cis' form to the extended 'trans' form."

  explanation = (
      "To be functionally accurate, the azobenzene-shaped picnic table must mimic the molecule's photoisomerization (shape-changing) properties:\n\n"
      "1. Azobenzene has a stable, extended form called the 'trans' isomer and a less stable, bent form called the 'cis' isomer.\n\n"
      f"2. {sunrise_action}\n\n"
      f"3. {sunset_action}\n\n"
      "This cycle would repeat daily, making the table a functional representation of a molecular machine."
  )
  print(explanation)

explain_azobenzene_table()