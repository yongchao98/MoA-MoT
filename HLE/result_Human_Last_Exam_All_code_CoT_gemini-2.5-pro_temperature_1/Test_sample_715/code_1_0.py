import pandas as pd
import requests
from bs4 import BeautifulSoup
import re

def clean_player_name(name):
    """
    Normalizes player names by removing non-alphanumeric characters,
    special symbols used by Baseball-Reference, and extra whitespace.
    """
    # Using unidecode would be ideal but it's not a standard library.
    # This simplified cleaning should work for the names in this specific problem.
    # B-R uses special characters for HOF, etc.
    name = re.sub(r'#|\*', '', name)
    # Normalize unicode characters like accents if possible without extra libraries
    # For this problem, direct matching after stripping suffices.
    return name.strip()

def get_players_from_game(url, team_name="Blue Jays"):
    """
    Parses a Baseball-Reference box score URL and returns a set of players for the specified team.
    """
    try:
        response = requests.get(url)
        response.raise_for_status()
    except requests.exceptions.RequestException as e:
        print(f"Error fetching {url}: {e}")
        return set()

    # Baseball-Reference wraps some tables in comments, need to remove them
    soup = BeautifulSoup(response.content, 'html.parser')
    comments = soup.find_all(string=lambda text: isinstance(text, BeautifulSoup.Comment))
    for comment in comments:
        comment_soup = BeautifulSoup(comment, 'html.parser')
        if comment_soup.find('table'):
            comment.replace_with(comment_soup)
    
    all_players = set()
    
    # Find tables that belong to the Blue Jays
    tables = soup.find_all('table', class_='sortable')
    for table in tables:
        caption = table.find('caption')
        if caption and team_name in caption.text:
            df = pd.read_html(str(table), header=0)[0]
            # Drop the last row which is 'Team Totals'
            df = df.iloc[:-1]
            # Player names are in the first column, labeled 'Batting' or 'Pitching'
            player_col_name = df.columns[0]
            # Add players to the set
            players = df[player_col_name].apply(clean_player_name)
            all_players.update(players)
            
    return all_players

def get_il_players(url):
    """
    Parses a Baseball-Reference transaction page to find players placed on the Injured List.
    """
    try:
        response = requests.get(url)
        response.raise_for_status()
    except requests.exceptions.RequestException as e:
        print(f"Error fetching {url}: {e}")
        return set()

    soup = BeautifulSoup(response.content, 'html.parser')
    transactions = soup.select('ul.page_index > li')
    
    il_players = set()
    # Regex to find names of players placed on the IL
    il_pattern = re.compile(r'placed\s(.*?)\son\s.*?\sinjured\slist', re.IGNORECASE)

    for trans in transactions:
        match = il_pattern.search(trans.text)
        if match:
            player_name = match.group(1).strip()
            il_players.add(clean_player_name(player_name))
            
    return il_players

def get_player_war(url):
    """
    Fetches and combines batting and pitching stats from a team's B-R page,
    returning a DataFrame with player names and their WAR.
    """
    try:
        response = requests.get(url)
        response.raise_for_status()
    except requests.exceptions.RequestException as e:
        print(f"Error fetching {url}: {e}")
        return pd.DataFrame()

    # Again, handle commented out tables
    html_content = response.text.replace('<!--', '').replace('-->', '')
    
    # Read batting and pitching tables by their IDs
    batting_df = pd.read_html(html_content, attrs={'id': 'team_batting'})[0]
    pitching_df = pd.read_html(html_content, attrs={'id': 'team_pitching'})[0]
    
    # Clean up and combine tables
    # Select only Name and WAR columns
    batting_war = batting_df[['Name', 'WAR']].copy()
    pitching_war = pitching_df[['Name', 'WAR']].copy()
    
    # Combine into a single DataFrame
    full_team_war = pd.concat([batting_war, pitching_war], ignore_index=True)
    
    # Clean player names
    full_team_war['Name'] = full_team_war['Name'].apply(clean_player_name)
    
    # Remove rows where 'Name' is 'Team Totals' or other non-player rows
    full_team_war = full_team_war[~full_team_war['Name'].str.contains('Total|Rank', case=False)]
    
    # Convert WAR to numeric, coercing errors
    full_team_war['WAR'] = pd.to_numeric(full_team_war['WAR'], errors='coerce')
    full_team_war.dropna(subset=['Name', 'WAR'], inplace=True)
    
    return full_team_war

# URLs for the required data
URL_GAME_1 = "https://www.baseball-reference.com/boxes/TOR/TOR201903280.shtml"
URL_GAME_FINAL = "https://www.baseball-reference.com/boxes/TOR/TOR201909290.shtml"
URL_TRANSACTIONS = "https://www.baseball-reference.com/teams/TOR/2019-transactions.shtml"
URL_TEAM_STATS = "https://www.baseball-reference.com/teams/TOR/2019.shtml"

# Step 1 & 2: Get players from the first and last games
players_game1 = get_players_from_game(URL_GAME_1)
players_game_final = get_players_from_game(URL_GAME_FINAL)

# Step 3: Get players who were on the Injured List
players_on_il = get_il_players(URL_TRANSACTIONS)

# Step 4: Find the intersection of all three sets to get candidate players
candidate_players = players_game1.intersection(players_game_final).intersection(players_on_il)

# Step 5: Get WAR for all players on the 2019 team
team_war_stats = get_player_war(URL_TEAM_STATS)

# Step 6: Filter the WAR stats for our candidate players
candidate_stats = team_war_stats[team_war_stats['Name'].isin(candidate_players)].copy()

# Step 7: Find the candidate with the highest WAR
if not candidate_stats.empty:
    winner = candidate_stats.loc[candidate_stats['WAR'].idxmax()]
    print(f"Player meeting all criteria: {winner['Name']}")
    print(f"2019 Baseball-Reference WAR: {winner['WAR']}")
else:
    print("No player found who meets all three requirements.")

final_answer_player = candidate_stats.loc[candidate_stats['WAR'].idxmax()]
print(f"<<<{final_answer_player['Name']}>>>")