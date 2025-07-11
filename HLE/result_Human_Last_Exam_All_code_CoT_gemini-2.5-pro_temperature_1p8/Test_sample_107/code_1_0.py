import ipaddress

def find_summary_acl():
    """
    Finds the smallest summary ACL entry for a given list of networks.
    """
    # Input networks in CIDR format
    network_strings = [
        "172.20.96.0/19",
        "172.20.128.0/19"
    ]

    # Convert strings to network objects
    try:
        networks = [ipaddress.ip_network(net_str) for net_str in network_strings]
    except ValueError as e:
        print(f"Error parsing network strings: {e}")
        return

    # Use collapse_addresses to find the minimal set of networks covering the space.
    # For this problem, it will be a single summary network.
    summary_networks = list(ipaddress.collapse_addresses(networks))

    if not summary_networks:
        print("Could not determine a summary network.")
        return
        
    # We expect a single summary network for this problem
    summary_net = summary_networks[0]
    
    # The network address of the summary route
    network_address = summary_net.network_address
    
    # The wildcard mask (hostmask) for the summary route
    wildcard_mask = summary_net.hostmask
    
    # Print the final ACL entry. The components of the entry are the network
    # address and the wildcard mask, fulfilling the requirement to output each number.
    print(f"{network_address} {wildcard_mask}")

if __name__ == "__main__":
    find_summary_acl()