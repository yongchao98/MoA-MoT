import ipaddress

def find_supernet_acl():
    """
    Calculates the smallest supernet (summary network) for a given list of networks
    and prints it as an access control list (ACL) entry with a wildcard mask.
    """
    # The networks provided by the user
    networks_str = ['172.20.96.0/19', '172.20.128.0/19']

    # Convert string representations to IPv4Network objects
    try:
        ip_networks = [ipaddress.ip_network(net) for net in networks_str]
    except ValueError as e:
        print(f"Error parsing network strings: {e}")
        return

    # collapse_addresses finds the minimal set of networks that covers the given networks.
    # In this case, it will return a single, summarized supernet.
    summarized_networks = list(ipaddress.collapse_addresses(ip_networks))

    # For this problem, we expect a single supernet as the result.
    if len(summarized_networks) == 1:
        supernet = summarized_networks[0]
        
        # The network address of the supernet
        network_address = supernet.network_address
        
        # The hostmask is the wildcard mask for an ACL entry
        wildcard_mask = supernet.hostmask
        
        print(f"The smallest ACL entry to match the networks is:")
        print(f"{network_address} {wildcard_mask}")
    else:
        print("Could not summarize the given networks into a single supernet.")
        print("Resulting networks:")
        for net in summarized_networks:
            print(f"{net.network_address} {net.hostmask}")

# Execute the function
find_supernet_acl()