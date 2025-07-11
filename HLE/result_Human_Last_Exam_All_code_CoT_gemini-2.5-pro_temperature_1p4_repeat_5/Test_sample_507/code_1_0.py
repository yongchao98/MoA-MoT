def find_linguistic_coincidences():
    """
    This script presents a solution to the linguistic puzzle by comparing
    words from Vietnamese and Persian (Farsi). These two languages are from
    different cultures in Asia with no direct historical contact.

    - Vietnamese: An Austroasiatic language spoken in Southeast Asia.
    - Persian (Farsi): An Indo-European language spoken in West Asia (Iran).

    The similarities, especially for "mom" and "dad", are common across many
    unrelated languages. However, the phonetic resemblance for "broom" is
    a more striking coincidence.
    """
    print("Comparing similar words from Vietnamese and Persian (Farsi):")
    print("-" * 60)

    # Data: English word, Vietnamese word, and Persian (Farsi) word
    word_comparisons = [
        {"english": "Mom", "vietnamese": "mẹ", "persian": "maman (مامان)"},
        {"english": "Dad", "vietnamese": "ba", "persian": "bābā (بابا)"},
        {"english": "Broom", "vietnamese": "chổi", "persian": "jāru (جارو)"}
    ]

    # Print the table header
    print(f"{'English':<10} | {'Vietnamese':<15} | {'Persian (Farsi)':<20}")
    print("="*60)

    # Print each word comparison
    for item in word_comparisons:
        print(f"{item['english']:<10} | {item['vietnamese']:<15} | {item['persian']:<20}")

    print("-" * 60)
    print("\nExplanation of similarities:")
    print("1. Dad: The Vietnamese 'ba' is nearly identical to the Persian 'bābā'.")
    print("2. Mom: The Vietnamese 'mẹ' and Persian 'maman' both start with the 'm' sound, a common feature of 'mother' words globally.")
    print("3. Broom: The Vietnamese 'chổi' and Persian 'jāru' show phonetic similarity. The initial sounds 'ch' [tɕ] and 'j' [dʒ] are very close, and both words feature a prominent vowel sound.")

find_linguistic_coincidences()
