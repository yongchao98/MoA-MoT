import collections

def analyze_emitter_lifetimes():
    """
    Analyzes iridium complexes to predict which will have shorter lifetimes in LECs.

    The stability of iridium emitters, which correlates with device lifetime,
    is highly sensitive to the substitution pattern on the cyclometalating ligands.
    
    A key principle is that substituents at the 'ortho' position (position 2)
    of the phenyl ring (relative to the Ir-C bond) introduce steric strain.
    This strain can destabilize the complex and weaken the Ir-C bond, leading
    to faster degradation and shorter lifetimes.
    """
    
    # Define the complexes based on the presence of ortho-Fluorine
    # True means an ortho-Fluorine is present, which leads to shorter lifetime.
    complexes = {
        1: {'has_ortho_fluorine': False},
        2: {'has_ortho_fluorine': False},
        3: {'has_ortho_fluorine': True},
        4: {'has_ortho_fluorine': True}
    }

    shorter_lifetime_complexes = []
    print("Analysis of Iridium Complexes for LEC Lifetime:")
    print("-" * 50)
    print("Principle: Fluorine atoms at the ortho-position of the cyclometalating\n"
          "phenyl ring cause steric hindrance, which destabilizes the complex\n"
          "and leads to a shorter operational lifetime.\n")

    print("Checking each complex for ortho-fluorination:")
    for number, properties in complexes.items():
        if properties['has_ortho_fluorine']:
            print(f" - Complex {number}: Contains an ortho-fluorine. Expected to be unstable.")
            shorter_lifetime_complexes.append(number)
        else:
            print(f" - Complex {number}: Lacks an ortho-fluorine. Expected to be relatively stable.")
    
    shorter_lifetime_complexes.sort()
    
    print("\nConclusion:")
    print("The complexes expected to show shorter lifetimes are those with ortho-fluorines.")
    
    # Print the final list of numbers as requested
    final_answer_string = " and ".join(map(str, shorter_lifetime_complexes))
    print(f"The identified complexes are {final_answer_string}.")

analyze_emitter_lifetimes()