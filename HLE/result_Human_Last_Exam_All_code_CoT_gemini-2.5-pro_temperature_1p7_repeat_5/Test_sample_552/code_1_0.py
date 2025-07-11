import sys
import io

# Buffer to capture the output
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

def find_cited_scientist():
    """
    Identifies the scientist whose theories were foundational to Einstein's 1905 Brownian motion paper.

    While Einstein's 1905 paper on Brownian motion famously lacks a formal bibliography,
    its title and arguments are explicitly based on the "Molecular-Kinetic Theory of Heat".
    This theory was primarily developed and championed by a specific physicist.
    Invoking the theory is a direct acknowledgment of its creator's work.
    """
    
    # The scientist most famously associated with the Molecular-Kinetic Theory of Heat
    scientist_name = "Ludwig Boltzmann"
    
    # The theory cited by name in Einstein's paper title
    theory_name = "Molecular-Kinetic Theory of Heat"
    
    # Print the explanation and the final answer
    print(f"Albert Einstein's 1905 paper on Brownian motion is explicitly based on the '{theory_name}'.")
    print("This theory was pioneered by, and is fundamentally the work of:")
    print(f"\n{scientist_name}")

find_cited_scientist()

# Restore original stdout and get the captured output
sys.stdout = old_stdout
output = captured_output.getvalue()

# Print the captured output to the user
print(output)

# Final Answer Extraction
final_answer_name = "Ludwig Boltzmann"
<<<final_answer_name>>>