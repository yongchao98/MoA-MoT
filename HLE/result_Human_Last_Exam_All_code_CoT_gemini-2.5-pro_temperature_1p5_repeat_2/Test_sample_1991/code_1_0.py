def create_cancer_phrase():
    """
    This function crafts and formats a phrase based on several specific constraints:
    1. Content: Describes abnormal cell growth (cancer) in the Luys' body.
    2. Syllables: Exactly fourteen.
    3. Rhyme/Scansion: Inspired by the 14-syllable Disney song title "Supercalifragilisticexpialidocious".
    4. Formatting: Each word capitalized.
    """
    
    # The crafted phrase:
    # - "Neoplastic Growth" refers to cancerous tumors.
    # - "In The Luys Body" specifies the location.
    # - "Atrocious" describes the condition and rhymes with "-docious".
    # - Syllable count: Ne-o-plas-tic (4) + Growth (1) + In (1) + The (1) + Luys (1) + Bo-dy (2) + Is (1) + A-tro-cious (3) = 14.
    base_phrase = "Neoplastic Growth In The Luys Body Is Atrocious"
    
    # Capitalize each word in the phrase.
    formatted_phrase = base_phrase.title()
    
    # Print the final result.
    print(formatted_phrase)

create_cancer_phrase()