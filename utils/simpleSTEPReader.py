import re
import numpy as np

_numeric_const_pattern = r"""
                        [-+]? # optional sign
                        (?:
                        (?: \d* \. \d+ ) # .1 .12 .123 etc 9.1 etc 98.1 etc
                        |
                        (?: \d+ \.? ) # 1. 12. 123. etc 1 12 123 etc
                        )
                        # followed by optional exponent part if desired
                        (?: [Ee] [+-]? \d+ ) ?
                        """

RE_NUMBERS=re.compile(_numeric_const_pattern, re.VERBOSE)
INSTANCE_DEFINITION_RE = re.compile("#(\d+)[^\S\n]?=[^\S\n]?(.*?)\((.*)\)[^\S\n]?;[\\r]?$")
#POINT_DEFINITION_RE = re.compile( ".*Point.[0-9]*\',\.*" )
POINT_DEFINITION_RE = re.compile( "\'.*\',\.*" )

def extractPointsFromSTEPFile(filename,reverse=False,mirrorCoord=None):
	fp = open(filename)
	lineDict={}
	i=0
	while True:
		line = fp.readline()
		if not line:
			break
		while (line.find(';') == -1): #its a multiline
			line = line.replace("\n","").replace("\r","") + fp.readline()
		match_instance_definition = INSTANCE_DEFINITION_RE.search(line)  # id,name,attrs
		if match_instance_definition:
			instance_id, entity_name, entity_attrs = match_instance_definition.groups()
			if entity_name == 'CARTESIAN_POINT':
				l=POINT_DEFINITION_RE.findall(entity_attrs)
				if len(l) >0:
					lineFloatList=RE_NUMBERS.findall( entity_attrs.replace(l[0],''))
					lineDict[i]=lineFloatList;i+=1;

	fp.close()
	rows=len(lineDict);
	cols=len(lineDict[0]);
	res=np.zeros((rows,cols))
	for i in range(0,rows):
		for j in range(0,cols):
			res[i][j]=float(lineDict[i][j])

	if np.all(res[0]==0.):
		res=res[1::]
	if mirrorCoord !=None:
		res[:,mirrorCoord]*=-1.


	if reverse:

		return res[::-1]
		
	else:
		return res


