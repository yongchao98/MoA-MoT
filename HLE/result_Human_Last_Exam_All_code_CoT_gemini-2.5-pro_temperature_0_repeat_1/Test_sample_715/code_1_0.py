import requests
from bs4 import BeautifulSoup
import re

def get_soup_from_url(url):
    """
    Fetches content from a URL and returns a BeautifulSoup object.
    It specifically handles the commented-out HTML tables on Baseball-Reference.
    """
    try:
        response = requests.get(url, headers={'User-Agent': 'Mozilla/5.0'})
        response.raise_for_status()
        # Baseball-Reference wraps tables in comments, so we extract the commented part
        comment_search = re.search(r'<!--(.*?)-->', response.text, re.DOTALL)
        if comment_search:
            content = comment_search.group(1)
        else:
            content = response.text
        return BeautifulSoup(content, 'html.parser')
    except requests.exceptions.RequestException as e:
        print(f"Error fetching URL {url}: {e}")
        return None

def get_players_from_box_score(soup):
    """Extracts a set of player names from a Blue Jays box score soup object."""
    players = set()
    if not soup:
        return players
    
    for table_id in ['TORbatting', 'TORpitching']:
        table = soup.find('table', id=table_id)
        if table:
            for link in table.find_all('a', href=re.compile(r'/players/')):
                players.add(link.text)
    return players

def get_il_players(soup):
    """Extracts a set of player names who were placed on the Injured List."""
    il_players = set()
    if not soup:
        return il_players
        
    for p_tag in soup.find_all('p'):
        text = p_tag.get_text(strip=True)
        if 'placed' in text.lower() and 'injured list' in text.lower():
            player_link = p_tag.find('a', href=re.compile(r'/players/'))
            if player_link:
                il_players.add(player_link.text)
    return il_players

def get_player_wars(soup):
    """Extracts a dictionary of player names and their WAR from a team stats soup object."""
    player_wars = {}
    if not soup:
        return player_wars

    for table_id in ['team_batting', 'team_pitching']:
        table = soup.find('table', id=table_id)
        if table:
            for row in table.find('tbody').find_all('tr'):
                player_cell = row.find('td', {'data-stat': 'player'})
                war_cell = row.find('td', {'data-stat': 'WAR'})
                
                if player_cell and war_cell and player_cell.a:
                    player_name = player_cell.a.text
                    try:
                        war_value = float(war_cell.text)
                        player_wars[player_name] = player_wars.get(player_name, 0) + war_value
                    except (ValueError, TypeError):
                        continue
    return player_wars

def solve_task():
    """Main function to execute the plan and find the player."""
    # URLs for the required Baseball-Reference pages
    URL_OPENING_DAY = 'https://www.baseball-reference.com/boxes/TOR/TOR201903280.shtml'
    URL_FINAL_DAY = 'https://www.baseball-reference.com/boxes/TOR/TOR201909290.shtml'
    URL_TRANSACTIONS = 'https://www.baseball-reference.com/teams/TOR/2019-transactions.shtml'
    URL_TEAM_STATS = 'https://www.baseball-reference.com/teams/TOR/2019.shtml'

    # Step 1: Get the three lists of players
    opening_day_soup = get_soup_from_url(URL_OPENING_DAY)
    opening_day_players = get_players_from_box_score(opening_day_soup)

    final_day_soup = get_soup_from_url(URL_FINAL_DAY)
    final_day_players = get_players_from_box_score(final_day_soup)

    transactions_soup = get_soup_from_url(URL_TRANSACTIONS)
    il_players = get_il_players(transactions_soup)

    # Step 2: Find players who meet all three criteria
    candidate_players = opening_day_players.intersection(final_day_players).intersection(il_players)

    if not candidate_players:
        print("No players were found who meet all three criteria.")
        return

    # Step 3: Get WAR for all players
    team_stats_soup = get_soup_from_url(URL_TEAM_STATS)
    player_wars = get_player_wars(team_stats_soup)

    # Step 4: Find the candidate with the highest WAR
    best_player = None
    highest_war = -float('inf')

    for player in candidate_players:
        war = player_wars.get(player)
        if war is not None and war > highest_war:
            highest_war = war
            best_player = player

    if best_player:
        print(f"The player who meets all three requirements and has the highest WAR is {best_player}.")
        print(f"This player's WAR in 2019 was: {highest_war}")
        print(f"<<<{best_player}>>>")
    else:
        print("Could not determine the player with the highest WAR among the candidates.")

if __name__ == '__main__':
    solve_task()