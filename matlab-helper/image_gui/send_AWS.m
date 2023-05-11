function send_AWS(datafolder,cycle_name,bucket)
cmd1=append('aws s3 cp e:\',datafolder,'\',cycle_name,' s3://',bucket,'/',datafolder,'/',cycle_name,' --storage-class DEEP_ARCHIVE');       
cmd2=append('aws s3 cp e:\',datafolder,'\','focus',cycle_name,' s3://',bucket,'/',datafolder,'/','focus',cycle_name,'  --storage-class DEEP_ARCHIVE');
cmd3=append('aws s3 cp e:\',datafolder,'\','dicfocus',cycle_name,' s3://',bucket,'/',datafolder,'/','dicfocus',cycle_name,' --recursive --storage-class DEEP_ARCHIVE');
cmd4=append('aws s3 cp e:\',datafolder,'\',cycle_name,'.csv',' s3://',bucket,'/',datafolder,'/','  --storage-class DEEP_ARCHIVE');
cmd5=append('aws s3 cp e:\',datafolder,'\','offset',cycle_name,'.csv',' s3://',bucket,'/',datafolder,'/',' --storage-class DEEP_ARCHIVE');
cmd6=append('aws s3 cp e:\',datafolder,'\','regoffset',cycle_name,'.csv',' s3://',bucket,'/',datafolder,'/','  --storage-class DEEP_ARCHIVE');
cmd7=append('aws s3 cp e:\',datafolder,'\','tiledregoffset',cycle_name,'.csv',' s3://',bucket,'/',datafolder,'/','  --storage-class DEEP_ARCHIVE');
system(cmd1);
system(cmd2);
system(cmd3)
system(cmd4);
system(cmd5);
system(cmd6)
system(cmd7)