def find_languages():
    """
    Identifies two languages based on a set of linguistic clues
    by filtering a pre-defined knowledge base.
    """

    # A knowledge base representing language properties.
    # This data is curated to demonstrate the logical deduction process.
    language_data = [
        {
            'name': 'Scottish Gaelic',
            'is_minority': True,
            'uses_k': False,
            'uses_w': False,
            'uses_à': True,
            'common_clusters': ['chd', 'rdb']
        },
        {
            'name': 'Northern Sámi',
            'is_minority': True,
            'uses_k': True,
            'uses_w': False,
            'uses_à': True,
            'common_clusters': ['ggj', 'skt', 'dj']
        },
        {
            'name': 'Italian',
            'is_minority': True, # e.g., in Slovenia
            'uses_k': True, # In loanwords
            'uses_w': True, # In loanwords
            'uses_à': True,
            'common_clusters': ['gli', 'gn', 'zz']
        },
        {
            'name': 'Welsh',
            'is_minority': True,
            'uses_k': False, # 'c' is used for the /k/ sound
            'uses_w': True,  # 'w' is a vowel
            'uses_à': False,
            'common_clusters': ['dd', 'll', 'ngh']
        }
    ]

    language_a = "Not Found"
    language_b = "Not Found"

    # Find Language 'a': no 'k' or 'w', has 'à', is a minority language.
    for lang in language_data:
        if (not lang['uses_k'] and
            not lang['uses_w'] and
            lang['uses_à'] and
            lang['is_minority']):
            language_a = lang['name']
            break

    # Find Language 'b': uses "ggj" and "skt", is a minority language.
    for lang in language_data:
        if ('ggj' in lang['common_clusters'] and
            'skt' in lang['common_clusters'] and
            lang['is_minority']):
            language_b = lang['name']
            break

    print(f"Language a: Based on the clues (no 'k' or 'w', uses 'à'), the language is {language_a}.")
    print(f"Language b: Based on the clues (uses 'ggj' and 'skt'), the language is {language_b}.")

# Run the function to solve the puzzle.
find_languages()