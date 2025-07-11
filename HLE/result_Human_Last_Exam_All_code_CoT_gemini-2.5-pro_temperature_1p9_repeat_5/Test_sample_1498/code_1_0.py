# The Creed of the Laymen ('Aqeedat al-'Awaam) is a 57-bayt (line) poem.
# This script identifies the specific range of abyat (lines) that list the names
# of the 25 prophets one is required to know.

# According to the standard numbering of the poem, the introduction to the list
# of prophets is in bayt 24. The names themselves begin on the next line.
start_bayt = 25

# The list of names is completed over the course of four abyat (lines).
# Bayt 25: Adam, Idris, Nuh, Hud, Salih, Ibrahim
# Bayt 26: Lut, Isma'il, Ishaq, Ya'qub, Yusuf, Ayyub
# Bayt 27: Shu'ayb, Harun, Musa, al-Yasa', Dhu al-Kifl, Dawud, Sulayman
# Bayt 28: Ilyas, Yunus, Zakariyya, Yahya, 'Isa, Taha (Muhammad)
end_bayt = 28

# The final print statement will clearly state the range as requested.
print(f"In the poem 'Aqeedat al-'Awaam', the names of the prophets are mentioned from bayt {start_bayt} to bayt {end_bayt}.")
