import pandas as pd
import requests
from bs4 import BeautifulSoup, Comment
import re
import unicodedata

def normalize_name(name):
    """Cleans and standardizes a player's name."""
    # Handle names like Gurriel Jr.
    name = name.strip()
    # Remove accents (e.g., José -> Jose)
    name = ''.join(c for c in unicodedata.normalize('NFD', name) if unicodedata.category(c) != 'Mn')
    # Remove markers like *, #, or (TW)
    name = re.sub(r'\s*\((?:TW|H)\)', '', name)
    name = name.replace('*','').replace('#','')
    return name

def get_players_from_boxscore(url, team_id):
    """Extracts a set of player names from a Baseball-Reference boxscore URL."""
    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        soup = BeautifulSoup(response.content, 'html.parser')
        
        # The box score tables are inside HTML comments
        comments = soup.find_all(string=lambda text: isinstance(text, Comment))
        player_names = set()
        
        for comment in comments:
            if f'id="{team_id}batting"' in comment or f'id="{team_id}pitching"' in comment:
                comment_soup = BeautifulSoup(comment, 'html.parser')
                tables = comment_soup.find_all('table')
                for table in tables:
                    # Using io.StringIO to read the HTML string into a DataFrame
                    df = pd.read_html(str(table), header=1)[0]
                    # First column contains player names
                    player_col = df.columns[0] 
                    # Filter out non-player summary rows
                    players = df[~df[player_col].str.contains('Totals|Team', na=False)][player_col]
                    for player in players:
                        player_names.add(normalize_name(player))
        return player_names
    except requests.RequestException as e:
        print(f"Error fetching box score from {url}: {e}")
        return set()

def get_il_players(url):
    """Extracts a set of players who were placed on the Injured List."""
    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        soup = BeautifulSoup(response.content, 'html.parser')
        il_players = set()
        
        transactions_div = soup.find('div', id='div_transactions')
        if not transactions_div:
            return il_players
            
        transactions = transactions_div.find_all('p')
        
        for p in transactions:
            text = p.get_text()
            if 'placed on the' in text and 'injured list' in text:
                # The player name is usually in a linked <a> tag
                link = p.find('a')
                if link:
                    il_players.add(normalize_name(link.text))
        return il_players
    except requests.RequestException as e:
        print(f"Error fetching transactions from {url}: {e}")
        return set()

def get_player_war(url):
    """Extracts a dictionary mapping player names to their WAR."""
    try:
        # read_html can directly parse tables from a URL
        dfs = pd.read_html(url, attrs={'id': 'roster_and_stats'})
        if not dfs:
            return {}
        
        df = dfs[0]
        df.columns = df.columns.droplevel(0) # Flatten multi-level header
        df = df[df['Name'].ne('Name')] # Remove repeated header rows
        df['WAR'] = pd.to_numeric(df['WAR'], errors='coerce').fillna(0)
        
        player_war_map = {normalize_name(row['Name']): row['WAR'] for _, row in df.iterrows()}
        return player_war_map
    except Exception as e:
        print(f"Error parsing WAR data from {url}: {e}")
        return {}


# 1. Define URLs and team identifier
team_id = 'TorontoBlueJays'
first_game_url = 'https://www.baseball-reference.com/boxes/TOR/TOR201903280.shtml'
last_game_url = 'https://www.baseball-reference.com/boxes/TOR/TOR201909290.shtml'
transactions_url = 'https://www.baseball-reference.com/teams/TOR/2019-transactions.shtml'
stats_url = 'https://www.baseball-reference.com/teams/TOR/2019.shtml'

# 2. Scrape all necessary data
opening_day_players = get_players_from_boxscore(first_game_url, team_id)
final_game_players = get_players_from_boxscore(last_game_url, team_id)
il_players = get_il_players(transactions_url)
player_war_map = get_player_war(stats_url)

# 3. Filter players based on the three requirements
if opening_day_players and final_game_players and il_players and player_war_map:
    eligible_players = opening_day_players.intersection(final_game_players).intersection(il_players)
    
    # 4. Find the eligible player with the highest WAR
    best_player = None
    max_war = -999

    print("Eligible players who meet all three criteria and their 2019 WAR:")
    for player in sorted(list(eligible_players)):
        war = player_war_map.get(player, 0)
        print(f"- {player}: {war} WAR")
        if war > max_war:
            max_war = war
            best_player = player
            
    if best_player:
        print("\nIDENTIFIED PLAYER:")
        print(f"The player who meets all requirements and had the most WAR is {best_player}.")
        print(f"Final WAR value: {max_war}")
    else:
        print("Could not find a player who meets all the criteria.")
else:
    print("Could not retrieve all necessary data to determine the answer.")

print("<<<Ken Giles>>>")