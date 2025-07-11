# First, ensure you have the geopy library installed:
# pip install geopy

from geopy.geocoders import Nominatim

def find_school_building_shape():
    """
    Finds the coordinates of the Pisa highschools and generates a
    satellite map URL to observe the shape of their buildings.
    """
    try:
        # Initialize the geolocator
        geolocator = Nominatim(user_agent="school_shape_finder_app")

        # The schools are adjacent on Via Benedetto Croce, Pisa.
        # We can use one address to center the map.
        address = "Liceo Scientifico Ulisse Dini, Pisa, Italy"
        
        print(f"Searching for coordinates of: {address}")
        location = geolocator.geocode(address, timeout=10)

        if location:
            # Extract latitude and longitude
            lat = location.latitude
            lon = location.longitude
            zoom = 19  # A close zoom level to see building details
            
            # Note: The numbers in the URL are the latitude, longitude, and zoom level.
            # Printing each number as requested by the prompt.
            print(f"Found Latitude: {lat}")
            print(f"Found Longitude: {lon}")
            print(f"Using Zoom Level: {zoom}")
            
            # Construct the Google Maps URL for satellite view
            # The format is https://www.google.com/maps/@{latitude},{longitude},{zoom}z/data=!3m1!1e3
            # The final '!3m1!1e3' argument forces satellite imagery.
            map_url = f"https://www.google.com/maps/@{lat},{lon},{zoom}z/data=!3m1!1e3"

            print("\n---------------------------------------------------------")
            print("By observing the satellite view at the URL below,")
            print("you can see the letter formed by the two school buildings.")
            print("URL for Satellite View:")
            print(map_url)
            print("---------------------------------------------------------")

        else:
            print("Could not find the location. Please check the address or your internet connection.")

    except ImportError:
        print("The 'geopy' library is required. Please install it by running: pip install geopy")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == '__main__':
    find_school_building_shape()