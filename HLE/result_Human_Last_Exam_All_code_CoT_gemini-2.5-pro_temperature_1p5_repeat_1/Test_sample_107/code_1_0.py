import ipaddress

def find_summary_acl_entry():
    """
    Calculates the smallest summary ACL entry for a list of networks.
    """
    # Define the input networks in CIDR notation
    network_strings = [
        '172.20.96.0/19',
        '172.20.128.0/19'
    ]

    try:
        # Convert the string representations into ipaddress.ip_network objects
        networks = [ipaddress.ip_network(net) for net in network_strings]

        # Use collapse_addresses to find the smallest possible summary route(s).
        # This function returns a generator, so we convert it to a list.
        summary_routes = list(ipaddress.collapse_addresses(networks))

        # For this specific case, we expect a single summary route.
        if summary_routes:
            print("The smallest summary ACL entry that covers the given networks is:")
            # An ACL entry consists of a network address and a wildcard mask.
            # The 'hostmask' attribute of an ip_network object is the wildcard mask.
            for route in summary_routes:
                network_address = route.network_address
                wildcard_mask = route.hostmask
                # The final result is the network address followed by the wildcard mask
                print(f"{network_address} {wildcard_mask}")
        else:
            print("Could not calculate a summary route for the given networks.")

    except ValueError as e:
        print(f"An error occurred: {e}. Please check the IP network formats.")

# Execute the function
find_summary_acl_entry()