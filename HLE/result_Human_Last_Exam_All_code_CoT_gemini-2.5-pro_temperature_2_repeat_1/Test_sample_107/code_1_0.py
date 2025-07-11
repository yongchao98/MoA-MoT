import ipaddress

# Define the list of IP networks to be summarized.
networks_to_summarize = [
    '172.20.96.0/19',
    '172.20.128.0/19'
]

try:
    # Use a list comprehension to convert the string representations 
    # into ipaddress.ip_network objects.
    ip_networks = [ipaddress.ip_network(net) for net in networks_to_summarize]

    # The ipaddress.collapse_addresses() function takes an iterable of networks
    # and returns the minimal set of networks that cover the same address space.
    # Since the input networks are contiguous, this will result in a single supernet.
    supernets = ipaddress.collapse_addresses(ip_networks)

    print("The smallest appropriate IP access control list entry is:")
    
    # Iterate through the resulting supernet(s) (in this case, just one).
    for net in supernets:
        # The network address is the base address of the network.
        network_address = net.network_address
        
        # The wildcard mask, used in ACLs, is the inverse of the subnet mask.
        # The ipaddress library provides this directly as the 'hostmask' attribute.
        wildcard_mask = net.hostmask
        
        # Print the final result in the format: <network address> <wildcard mask>
        print(f"{network_address} {wildcard_mask}")

except ValueError as e:
    print(f"An error occurred during IP address processing: {e}")
except Exception as e:
    print(f"An unexpected error occurred: {e}")
