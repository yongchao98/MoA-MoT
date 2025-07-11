import requests
from bs4 import BeautifulSoup
from bs4 import Comment

def get_players_from_box_score(url):
    """
    Scrapes a Baseball-Reference box score URL to get a set of player names for a team.
    The required tables are often hidden inside HTML comments, so this function parses them.
    """
    players = set()
    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        soup = BeautifulSoup(response.content, 'html.parser')
        
        # Find all commented-out sections in the HTML
        comments = soup.find_all(string=lambda text: isinstance(text, Comment))
        
        for comment in comments:
            # Check if the comment contains the team's batting or pitching tables
            if 'id="Torontobp_batting"' in comment or 'id="Torontobp_pitching"' in comment:
                box_soup = BeautifulSoup(comment, 'html.parser')
                tables = box_soup.find_all('table')
                for table in tables:
                    for row in table.find('tbody').find_all('tr'):
                        # Skip header rows that are sometimes included in the tbody
                        if 'thead' not in row.get('class', []):
                            player_cell = row.find('th')
                            if player_cell and player_cell.find('a'):
                                players.add(player_cell.find('a').text)
    except requests.exceptions.RequestException as e:
        print(f"Error fetching box score at {url}: {e}")
    return players

def get_il_players(team_abbr, year):
    """
    Scrapes the team's transaction page for a given year to find all players
    placed on the Injured List.
    """
    url = f"https://www.baseball-reference.com/teams/{team_abbr}/{year}-transactions.shtml"
    il_players = set()
    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        soup = BeautifulSoup(response.content, 'html.parser')
        
        # Transactions are stored in <p> tags within the main content div
        transaction_paragraphs = soup.find('div', id='content').find_all('p')
        
        for p in transaction_paragraphs:
            # Look for the specific text indicating an IL placement
            if "placed on the" in p.text and "injured list" in p.text:
                player_link = p.find('a')
                if player_link:
                    il_players.add(player_link.text)
    except requests.exceptions.RequestException as e:
        print(f"Error fetching transactions at {url}: {e}")
    return il_players

def get_player_war(team_abbr, year):
    """
    Scrapes the main team page for a given year to get the WAR for each player.
    It combines batting and pitching WAR for players who may have done both.
    """
    url = f"https://www.baseball-reference.com/teams/{team_abbr}/{year}.shtml"
    player_war = {}
    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        soup = BeautifulSoup(response.content, 'html.parser')

        # --- Scrape Batting WAR (this table is not in a comment) ---
        batting_table = soup.find('table', id='team_batting')
        for row in batting_table.find('tbody').find_all('tr'):
            if 'thead' in row.get('class', []):
                continue
            name_cell = row.find('td', {'data-stat': 'player'})
            war_cell = row.find('td', {'data-stat': 'WAR'})
            if name_cell and war_cell and name_cell.find('a'):
                name = name_cell.find('a').text
                try:
                    war = float(war_cell.text)
                    player_war[name] = player_war.get(name, 0) + war
                except (ValueError, TypeError):
                    continue

        # --- Scrape Pitching WAR (this table is in a comment) ---
        comments = soup.find_all(string=lambda text: isinstance(text, Comment))
        for comment in comments:
            if 'id="team_pitching"' in str(comment):
                pitching_soup = BeautifulSoup(str(comment), 'html.parser')
                pitching_table = pitching_soup.find('table', id='team_pitching')
                if pitching_table:
                    for row in pitching_table.find('tbody').find_all('tr'):
                        if 'thead' in row.get('class', []):
                            continue
                        name_cell = row.find('td', {'data-stat': 'player'})
                        war_cell = row.find('td', {'data-stat': 'WAR'})
                        if name_cell and war_cell and name_cell.find('a'):
                            name = name_cell.find('a').text
                            try:
                                war = float(war_cell.text)
                                player_war[name] = player_war.get(name, 0) + war
                            except (ValueError, TypeError):
                                continue
    except requests.exceptions.RequestException as e:
        print(f"Error fetching player WAR at {url}: {e}")
    return player_war

if __name__ == "__main__":
    # 1. Get players from the first game of the 2019 season
    first_game_url = "https://www.baseball-reference.com/boxes/TOR/TOR201903280.shtml"
    first_game_players = get_players_from_box_score(first_game_url)

    # 2. Get players from the last game of the 2019 season
    last_game_url = "https://www.baseball-reference.com/boxes/TOR/TOR201909290.shtml"
    last_game_players = get_players_from_box_score(last_game_url)

    # 3. Get all players who were placed on the IL in 2019
    il_players = get_il_players('TOR', 2019)

    # 4. Find the players who meet all three criteria
    eligible_players = first_game_players.intersection(last_game_players).intersection(il_players)
    
    if not eligible_players:
        print("No players found who meet all three criteria.")
    else:
        # 5. Get the WAR for all players on the 2019 team
        all_player_war = get_player_war('TOR', 2019)
        
        highest_war = -999.0
        top_player_name = None

        # 6. Iterate through eligible players to find the one with the highest WAR
        print("Eligible players who meet all three criteria:")
        for player in sorted(list(eligible_players)):
            war = all_player_war.get(player, 0)
            print(f"- {player} (WAR: {war})")
            if war > highest_war:
                highest_war = war
                top_player_name = player

        print("\n" + "="*40)
        print("Final Answer:")
        if top_player_name:
            print(f"The player with the highest WAR is {top_player_name} with a WAR of {highest_war}.")
        else:
            print("Could not determine the player with the highest WAR.")
        print("="*40)