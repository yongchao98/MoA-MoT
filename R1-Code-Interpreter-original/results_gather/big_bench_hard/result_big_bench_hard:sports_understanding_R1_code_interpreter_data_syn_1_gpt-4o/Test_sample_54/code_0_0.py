import requests
from bs4 import BeautifulSoup

def search_ryan_tannehill_baseball():
    query = "Ryan Tannehill baseball career"
    url = f"https://www.google.com/search?q={query}"
    headers = {
        "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/58.0.3029.110 Safari/537.3"
    }
    response = requests.get(url, headers=headers)
    soup = BeautifulSoup(response.text, "html.parser")
    results = soup.find_all('div', class_='BNeawe vvjwJb AP7Wnd')
    
    for result in results[:5]:  # Check the first 5 results
        print(result.get_text())

search_ryan_tannehill_baseball()