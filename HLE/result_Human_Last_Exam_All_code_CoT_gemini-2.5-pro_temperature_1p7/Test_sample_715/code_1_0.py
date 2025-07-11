import pandas as pd
import requests
from bs4 import BeautifulSoup, Comment
import re
import unicodedata
import io

def normalize_name(name):
    """
    Normalizes player names for consistent matching by handling unicode,
    suffixes (like Jr.), and periods.
    """
    if not isinstance(name, str):
        return ""
    # Handle unicode characters like accents
    try:
        nfkd_form = unicodedata.normalize('NFKD', name)
        name = u"".join([c for c in nfkd_form if not unicodedata.combining(c)])
    except TypeError:
        pass # Handle potential errors with non-string inputs gracefully
    # Remove suffixes like Jr., Sr., II, III, IV
    name = re.sub(r'\s+(Jr|Sr|II|III|IV)\.?$', '', name, flags=re.IGNORECASE)
    # Remove periods (e.g., from initials like A.J.)
    name = name.replace('.', '')
    return name.strip()

def get_players_from_game(url):
    """
    Scrapes a Baseball-Reference box score URL and returns a set of player names for the Blue Jays.
    """
    players = set()
    try:
        response = requests.get(url, timeout=15)
        response.raise_for_status()
        soup = BeautifulSoup(response.content, 'html.parser')
        
        # The team tables are easily identified by an ID containing 'TorontoBlueJays'
        team_tables = soup.find_all('table', id=re.compile(r'TorontoBlueJays(batting|pitching)'))
        for table in team_tables:
            for a in table.find_all('a', href=re.compile(r'/players/')):
                players.add(normalize_name(a.text))
    except requests.exceptions.RequestException as e:
        print(f"Error fetching game data from {url}: {e}")
    return players

def get_il_players(url):
    """
    Scrapes the team transaction page to find players placed on the Injured List.
    """
    il_players = set()
    try:
        response = requests.get(url, timeout=15)
        response.raise_for_status()
        soup = BeautifulSoup(response.content, 'html.parser')
        
        transactions_div = soup.find('div', id='all_transactions')
        if transactions_div:
            transactions = transactions_div.find_all('p')
            for trans in transactions:
                # Look for the specific phrase "placed on ... injured list"
                if re.search(r'placed on .* injured list', trans.text, re.IGNORECASE):
                    player_link = trans.find('a', href=re.compile(r'/players/'))
                    if player_link:
                        il_players.add(normalize_name(player_link.text))
    except requests.exceptions.RequestException as e:
        print(f"Error fetching IL data from {url}: {e}")
    return il_players

def get_player_war(url):
    """
    Scrapes the main team stats page (which has its primary table in a comment)
    to get the WAR for each player. It returns a dictionary mapping normalized
    player names to their WAR.
    """
    player_war_map = {}
    try:
        response = requests.get(url, timeout=15)
        response.raise_for_status()
        soup = BeautifulSoup(response.content, 'html.parser')
        
        # The main data table is embedded within an HTML comment
        comment = soup.find(string=lambda text: isinstance(text, Comment) and "team_batting" in text)
        if comment:
            # Use pandas to parse the HTML table from the comment string
            df = pd.read_html(io.StringIO(comment), attrs={'id': 'team_batting'})[0]
            
            df = df[df['Rk'].ne('Rk') & df['Player'].notna()].copy()
            df = df.dropna(subset=['Rk'])
            
            df['Player_norm'] = df['Player'].apply(normalize_name)
            df['WAR_val'] = pd.to_numeric(df['WAR'], errors='coerce').fillna(0.0)
            
            player_war_map = df.set_index('Player_norm')['WAR_val'].to_dict()

    except Exception as e:
        print(f"An error occurred: {e}")
        
    return player_war_map

def solve_blue_jays_riddle():
    # URLs for the 2019 Toronto Blue Jays
    first_game_url = 'https://www.baseball-reference.com/boxes/TOR/TOR201903280.shtml'
    last_game_url = 'https://www.baseball-reference.com/boxes/TOR/TOR201909290.shtml'
    transactions_url = 'https://www.baseball-reference.com/teams/TOR/2019-transactions.shtml'
    stats_url = 'https://www.baseball-reference.com/teams/TOR/2019.shtml'

    # 1. Get player lists based on the criteria
    players_first_game = get_players_from_game(first_game_url)
    players_last_game = get_players_from_game(last_game_url)
    players_on_il = get_il_players(transactions_url)
    player_war_map = get_player_war(stats_url)

    # 2. Find the set of players who meet all three criteria
    qualified_players = players_first_game.intersection(players_last_game).intersection(players_on_il)

    # 3. Find the qualified player with the highest WAR
    if not qualified_players or not player_war_map:
        print("Could not find any players matching all criteria or failed to retrieve data.")
        return

    best_player_name = ""
    max_war = -999

    for player_name in qualified_players:
        if player_name in player_war_map:
            current_war = player_war_map[player_name]
            if current_war > max_war:
                max_war = current_war
                best_player_name = player_name

    # 4. Print the final answer
    if best_player_name:
        print("The player who meets all criteria and has the highest WAR is:")
        print(f"Player Name: {best_player_name}")
        print(f"2019 WAR: {max_war}")
    else:
        print("Could not determine the player with the highest WAR from the qualified candidates.")

if __name__ == '__main__':
    solve_blue_jays_riddle()
<<<Ken Giles>>>