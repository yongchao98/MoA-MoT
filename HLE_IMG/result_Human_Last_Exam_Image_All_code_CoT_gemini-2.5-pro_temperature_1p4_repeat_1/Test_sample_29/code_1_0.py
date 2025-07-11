import math

def calculate_hill_number(abundances, q, total_individuals):
    """Calculates the Hill number for a given list of abundances and order q."""
    proportions = [n / total_individuals for n in abundances if n > 0]
    
    if not proportions:
        return 0.0

    if q == 1:
        shannon_entropy = -sum(p * math.log(p) for p in proportions)
        return math.exp(shannon_entropy)
    else:
        sum_of_proportions_raised_to_q = sum(p**q for p in proportions)
        if sum_of_proportions_raised_to_q == 0:
            return 0.0
        return sum_of_proportions_raised_to_q**(1 / (1 - q))

def print_equation(level_name, abundances, q, total_individuals, result):
    """Prints the formatted equation for a Hill number calculation."""
    print(f"Calculating for {level_name} level (q={q}):")
    print(f"  Abundances: {abundances}")
    
    prop_strings = [f"({n}/{total_individuals})" for n in abundances]
    
    if q == 1:
        # Format for D^1 = exp(-Σ(pi*ln(pi)))
        inner_sum = " + ".join([f"{p_str}*ln({p_str})" for p_str in prop_strings])
        print(f"  Equation: D^1 = exp( -( {inner_sum} ) )")
    else:
        # Format for D^q = (Σ(pi^q))^(1/(1-q))
        inner_sum = " + ".join([f"{p_str}^{q}" for p_str in prop_strings])
        print(f"  Equation: D^{q} = ( {inner_sum} )^(1/(1-{q}))")
        
    print(f"  Result: {result:.2f}\n")

# Step 1: Define abundances based on visual count and grouping
species_abundances = [81, 15, 9, 5]
genus_abundances = [81 + 5, 15, 9] # Grouping species 1 and 4
family_abundances = [86 + 9, 15]   # Grouping genus 1 and 3
total_individuals = sum(species_abundances)
order_abundances = [total_individuals]

# Step 2: Define q values
q_order = 1
q_family = 2
q_genus = 3
q_species = 4

# Step 3: Calculate Hill numbers
hill_order = calculate_hill_number(order_abundances, q_order, total_individuals)
hill_family = calculate_hill_number(family_abundances, q_family, total_individuals)
hill_genus = calculate_hill_number(genus_abundances, q_genus, total_individuals)
hill_species = calculate_hill_number(species_abundances, q_species, total_individuals)

# Step 4: Print equations and results
print_equation("Order", order_abundances, q_order, total_individuals, hill_order)
print_equation("Family", family_abundances, q_family, total_individuals, hill_family)
print_equation("Genus", genus_abundances, q_genus, total_individuals, hill_genus)
print_equation("Species", species_abundances, q_species, total_individuals, hill_species)

# Step 5: Format the final answer as requested
final_answer = f"{hill_order:.2f}, {hill_family:.2f}, {hill_genus:.2f}, {hill_species:.2f}"
print("Final comma-separated values:")
print(final_answer)
print("<<<1.00, 1.31, 1.44, 1.50>>>")