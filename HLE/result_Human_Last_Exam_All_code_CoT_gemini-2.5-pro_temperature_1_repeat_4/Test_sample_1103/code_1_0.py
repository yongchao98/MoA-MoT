import requests
import re

def find_discriminant_count_for_class_number(class_number):
    """
    Finds the number of negative fundamental discriminants for a given class number
    by querying the LMFDB (L-Functions and Modular Forms Database).
    """
    # The URL for imaginary quadratic fields on LMFDB, filtered by class number.
    url = f"https://www.lmfdb.org/NumberField/?hst=1&class_number={class_number}"
    
    print(f"Querying database for class number = {class_number}...")
    
    try:
        # Use a user-agent to mimic a browser, which can help avoid access issues.
        headers = {
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/58.0.3029.110 Safari/537.3'
        }
        response = requests.get(url, headers=headers)
        # Raise an HTTPError for bad responses (4xx or 5xx)
        response.raise_for_status()
        
        html_content = response.text
        
        # The page displays the total count in a format like "Displaying 1 to 50 of 50".
        # We use a regular expression to find and extract this total count.
        match = re.search(r"Displaying [0-9,]+ to [0-9,]+ of ([0-9,]+)", html_content)
        
        if match:
            # The extracted number might have commas, so we remove them before converting to an integer.
            count_str = match.group(1).replace(',', '')
            count = int(count_str)
            print(f"The equation is: Number of discriminants for class number {class_number} = {count}")
        elif "No NumberField labels match your query" in html_content:
            # This handles cases where the class number has no associated discriminants.
            count = 0
            print(f"The equation is: Number of discriminants for class number {class_number} = {count}")
        else:
            print("Could not find the count. The website's layout may have changed.")

    except requests.exceptions.RequestException as e:
        print(f"An error occurred while trying to fetch data from LMFDB: {e}")

# For the Gauss class number problem, we want the count for class number 48.
find_discriminant_count_for_class_number(48)