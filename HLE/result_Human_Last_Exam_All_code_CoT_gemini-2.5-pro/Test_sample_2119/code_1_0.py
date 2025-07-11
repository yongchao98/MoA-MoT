def find_minister_from_historical_context():
    """
    This function analyzes historical context to identify the minister's name.
    """
    # The text describes the banning of a tragedy about "Alexei".
    # This refers to plays about Tsarevich Alexei Petrovich, son of Peter the Great.
    # One such famous play is Dmitry Merezhkovsky's "Peter and Alexis".
    # The minister-like figure described as "seeking to restore order" by banning art
    # aligns perfectly with the historical role and conservative ideology of
    # Konstantin Pobedonostsev, the Ober-Procurator of the Most Holy Synod.
    # His position gave him ministerial power over cultural and religious censorship.

    historical_record = {
        "play_subject": "Tsarevich Alexei",
        "key_censor_title": "Ober-Procurator of the Holy Synod",
        "censor_name": "Konstantin Pobedonostsev",
        "ideology": "conservative, reactionary, 'restoring order'"
    }

    minister_full_name = historical_record.get("censor_name")
    
    if minister_full_name:
        # The last name is the second part of the full name.
        minister_last_name = minister_full_name.split()[1]
        print(minister_last_name)
    else:
        print("Minister not found in records.")

find_minister_from_historical_context()