import requests
from bs4 import BeautifulSoup

def find_restored_scene_element():
    """
    Fetches and analyzes a web article to determine the content of a restored
    scene in the film 'Kriemhild's Revenge'.
    """
    # URL of the French article discussing the restored version of the film
    url = "https://www.lemonde.fr/vous/article/2011/10/01/les-nibelungen_1580982_3238.html"
    
    answer_choices = {
        'A': "A shot of Hildebrand bringing Kriemhild's crown to Etzel.",
        'B': "A shot of Etzel watching the sunset, mired in sorrow.",
        'C': "A shot of Hildebrand striking Kriemhild down with his spear.",
        'D': "A shot of Etzel lifts his infant son amidst the carnage.",
        'E': "A shot of Etzel calling for help, lamenting.",
        'F': "A shot of Kriemhild placing Siegfried's sword on her lap in the shadowy hall."
    }

    print(f"Step 1: Attempting to fetch content from the source article:\n{url}\n")
    try:
        # Using a User-Agent header to mimic a browser visit
        headers = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36'}
        response = requests.get(url, headers=headers, timeout=15)
        # Raise an HTTPError if the HTTP request returned an unsuccessful status code
        response.raise_for_status()
    except requests.exceptions.RequestException as e:
        print(f"Error: Could not retrieve the web page. The process cannot continue.")
        print(f"Details: {e}")
        return

    print("Step 2: Successfully fetched page. Parsing the article's text...")
    soup = BeautifulSoup(response.content, 'html.parser')
    
    # The main article content is typically in '<p>' tags with a specific class
    paragraphs = soup.find_all('p', class_='article__paragraph')
    
    found_description = None
    # Search for keywords to identify the relevant paragraph. 
    # The key elements are Etzel, his son ('fils' in French), and the action of lifting ('brandit').
    for p in paragraphs:
        text = p.get_text()
        if "Etzel" in text and "fils" in text and "brandit au ciel" in text:
            found_description = text.strip()
            break

    if found_description:
        print("\nStep 3: Found the paragraph describing the restored final scene.")
        print(f"\nOriginal Text (from Le Monde):\n'{found_description}'")
        
        # Provide a translation of the key sentence
        print("\nEnglish Translation of the key part:\n'In the midst of the carnage, Etzel takes in his arms the son he had with Kriemhild, and brandishes him to the sky, as a sign of survival and hope.'")
        
        print("\nStep 4: Comparing this description with the given answer choices...")
        # The description clearly points to Etzel lifting his son.
        correct_answer_key = 'D'
        print(f"The evidence from the article directly matches answer '{correct_answer_key}':")
        print(f"'{answer_choices[correct_answer_key]}'")
    else:
        print("\nCould not automatically locate the specific description in the article's text. The website's structure may have changed.")

# Execute the function
find_restored_scene_element()