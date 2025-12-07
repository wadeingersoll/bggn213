#Open an SSH client.

#Locate your private key file. The key used to launch this instance is bggn213_wade_key.pem

#Run this command, if necessary, to ensure your key is not publicly viewable.
#chmod 400 "bggn213_wade_key.pem"

#Connect to your instance using its Public DNS:
#  ec2-54-214-110-23.us-west-2.compute.amazonaws.com

#Example:
  
#  ssh -i ~/Downloads/bggn213_wade_key.pem ubuntu@ec2-54-214-110-23.us-west-2.compute.amazonaws.com


# download quant files:

# scp -r -i ~/Downloads/bggn213_wade_key.pem ubuntu@ec2-54-214-110-23.us-west-2.compute.amazonaws.com:~/*_quant .

