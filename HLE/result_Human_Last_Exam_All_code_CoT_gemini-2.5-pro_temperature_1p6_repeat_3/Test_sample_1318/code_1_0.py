# Experiment 1 Data: Plaque-Forming Units (cfu/ul)
cfu_data = {
    "no_RP": {
        "wt": 100000,
        "deltaXY": 100000,
    },
    "with_RP": {
        "wt": 80000,
        "deltaXY": 40000,
    }
}

# Experiment 2 Data: Detection of 500 Da molecule at 60 minutes
mass_spec_data = {
    "Sample1_RP_wt": True,      # with RP, with XY
    "Sample2_RP_deltaXY": False,  # with RP, without XY
    "Sample3_noRP_wt": False,     # without RP, with XY
    "Sample4_noRP_deltaXY": False # without RP, without XY
}

print("--- Analysis of Experimental Data ---")

# Analysis of Statement H:
# "System RP increases the resistance of the bacteria against phageDE3 because
# the enzymes XY1 or XY2 can synthesize their products only in the presence of the system RP."

print("\nStep 1: Evaluate the first part of statement H -> 'System RP increases the resistance of the bacteria against phageDE3'")

# We check this by comparing the success of a phage (e.g., deltaXY) with and without the RP system.
cfu_no_rp = cfu_data["no_RP"]["deltaXY"]
cfu_with_rp = cfu_data["with_RP"]["deltaXY"]
is_resistance_increased = cfu_with_rp < cfu_no_rp

print(f"Comparing phageDE3-deltaXY infection:")
print(f"CFU in bacteria without RP system: {cfu_no_rp}")
print(f"CFU in bacteria with RP system: {cfu_with_rp}")
print(f"Since {cfu_with_rp} < {cfu_no_rp}, the RP system reduces phage success.")
print(f"Conclusion: The statement that RP increases resistance is TRUE.")


print("\nStep 2: Evaluate the second part of statement H -> '...the enzymes XY1 or XY2 can synthesize their products only in the presence of the system RP'")

# We check if the product appears when XY is present, but only if RP is also present.
# We compare Sample 1 (with RP, with XY) vs Sample 3 (without RP, with XY)
product_with_rp_and_xy = mass_spec_data["Sample1_RP_wt"]
product_without_rp_and_xy = mass_spec_data["Sample3_noRP_wt"]

print("Comparing production of 500 Da molecule for phageDE3-wt (with XY enzymes):")
print(f"Detection in bacteria with RP system: {product_with_rp_and_xy}")
print(f"Detection in bacteria without RP system: {product_without_rp_and_xy}")
print("Conclusion: The product is synthesized only when both the XY enzymes and the RP system are present. This statement is TRUE.")

print("\n--- Final Conclusion ---")
print("Both parts of statement H are supported by the experimental data. It provides the most complete explanation of the phenomena observed.")
print("The final answer is H.")
<<<H>>>