import requests
import re

def find_human_let7_members():
    """
    Fetches data from miRBase and counts the number of human let-7 family members.
    """
    # 1. Fetch the latest miRNA data from miRBase.
    url = "https://www.mirbase.org/ftp/CURRENT/miRNA.dat"
    print(f"Fetching data from {url}...")

    try:
        response = requests.get(url)
        response.raise_for_status()  # Raise an exception for bad status codes
        data = response.text
        print("Data fetched successfully.\n")
    except requests.exceptions.RequestException as e:
        print(f"Error fetching data: {e}")
        # As a fallback, use the known number if the database is inaccessible.
        print("Could not fetch live data. Based on current biological consensus, there are 9 members.")
        print("These are: let-7a, let-7b, let-7c, let-7d, let-7e, let-7f, let-7g, let-7i, and miR-98.")
        print("\nTotal number of members: 9")
        return

    # Use a set to store unique let-7 family members to avoid duplicates.
    let7_family_members = set()

    # Each entry in the file is separated by "//"
    entries = data.strip().split('//\n')

    # 2. Parse the data to find all entries for Homo sapiens (human).
    for entry in entries:
        if 'OS   Homo sapiens' in entry:
            # Find the ID line, e.g., "ID   hsa-let-7a-1   standard; RNA; H sapiens; 83 BP."
            id_match = re.search(r'^ID\s+([^\s]+)', entry, re.MULTILINE)
            if id_match:
                mirna_id = id_match.group(1)

                # 3. Identify all unique members of the let-7 family.
                # This includes 'hsa-let-7' followed by a letter and 'hsa-mir-98'.
                let7_match = re.match(r'hsa-let-7([a-z])', mirna_id)
                mir98_match = re.match(r'hsa-mir-98', mirna_id)

                if let7_match:
                    # Extract the core family member name, e.g., 'let-7a' from 'hsa-let-7a-1'
                    member_name = f"let-7{let7_match.group(1)}"
                    let7_family_members.add(member_name)
                elif mir98_match:
                    # Add miR-98 to the family
                    let7_family_members.add("miR-98")

    if not let7_family_members:
        print("Could not find any human let-7 family members. The database format might have changed.")
        return

    # 4. Count and display the members.
    sorted_members = sorted(list(let7_family_members))
    count = len(sorted_members)

    print("The following members of the let-7 family have been identified in humans:")
    # Output each member that contributes to the final count
    for member in sorted_members:
        print(f"1 # ({member})")

    print(f"\nEquation: {' + '.join(['1'] * count)} = {count}")
    print(f"\nTo date, {count} members of the let-7 family have been identified in humans.")

# Execute the function
find_human_let7_members()