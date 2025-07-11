# The task is to find the earliest known date recorded by a pre-Columbian civilization
# in an aboriginal writing system. The answer is found on Stela 2 from Chiapa de Corzo.

# Define the components of the Long Count date.
baktun = 7
katun = 16
tun = 3
winal = 2
kin = 13

# Define other relevant information.
artifact_name = "Stela 2"
location = "Chiapa de Corzo, Mexico"
gregorian_equivalent = "December 7, 36 BC"

# Print the final answer, detailing the date and its components.
print(f"The earliest known date from a pre-Columbian civilization is found on {artifact_name} from {location}.")
print(f"The date is written in the Long Count system as:")
print(f"The full date is: {baktun}.{katun}.{tun}.{winal}.{kin}")
print("\nThis breaks down into the following units of time:")
print(f"- {baktun} B'ak'tuns")
print(f"- {katun} K'atuns")
print(f"- {tun} Tuns")
print(f"- {winal} Winals")
print(f"- {kin} K'ins")
print(f"\nThis date corresponds to {gregorian_equivalent} in the Gregorian calendar.")