def analyze_bansenshukai_pattern():
    """
    Analyzes the Kunoichi no Jutsu pattern from the Bansenshukai.

    The pattern uses '⬤' for known kanji and '○' for missing kanji.
    This script will count each type and present a summary equation.
    """
    pattern_str = "⬤○○⬤⬤⬤⬤○⬤⬤⬤⬤⬤"

    # Count the occurrences of each symbol
    known_kanji_count = pattern_str.count('⬤')
    missing_kanji_count = pattern_str.count('○')
    total_symbols = len(pattern_str)

    print("Analyzing the Kunoichi no Jutsu Pattern:")
    print(f"Pattern: {pattern_str}\n")

    print(f"Number of known/existing kanji (⬤): {known_kanji_count}")
    print(f"Number of unknown/missing kanji (○): {missing_kanji_count}")
    print(f"Total number of symbols: {total_symbols}\n")

    # Display the final count as a simple equation as requested
    print("Final Equation:")
    print(f"{known_kanji_count} (Known) + {missing_kanji_count} (Missing) = {total_symbols} (Total Symbols)")

if __name__ == "__main__":
    analyze_bansenshukai_pattern()