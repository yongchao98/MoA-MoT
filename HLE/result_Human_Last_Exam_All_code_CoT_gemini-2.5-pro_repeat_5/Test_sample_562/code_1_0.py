import requests
from bs4 import BeautifulSoup
import re

def find_radical_name():
    """
    Finds the bushu/radical for the character 謄 by scraping Jisho.org
    and prints its name in Hepburn transcription.
    """
    kanji = '謄'
    
    try:
        # The URL for the specific kanji page on Jisho.org
        url = f"https://jisho.org/search/{kanji}%20%23kanji"
        headers = {
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36'
        }

        # Send a request to the website
        response = requests.get(url, headers=headers)
        response.raise_for_status()

        # Parse the HTML content
        soup = BeautifulSoup(response.text, 'html.parser')

        # Find the 'dt' element with the text 'Radical:'
        # Using regex to ignore potential whitespace
        radical_label = soup.find('dt', string=re.compile(r'\s*Radical:\s*'))

        if radical_label:
            # The radical name is in the next 'dd' element
            radical_data = radical_label.find_next_sibling('dd')
            if radical_data:
                # The text is structured like: '肉 nikuniku, nikuzuki'
                # We extract the full text, strip whitespace
                full_text = radical_data.get_text(strip=True)
                
                # The first character is the radical symbol, we remove it
                readings_part = full_text[1:]
                
                # Split readings by comma and take the last one, which is often the most specific name
                radical_name = readings_part.split(',')[-1].strip()
                
                # Capitalize the first letter as in the example "Onnahen"
                final_answer = radical_name.capitalize()
                
                print(final_answer)
                return

        print("Could not find the radical information.")

    except requests.exceptions.RequestException as e:
        print(f"Error fetching data: {e}")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == '__main__':
    find_radical_name()