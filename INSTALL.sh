#!/usr/bin/env sh

# Check pipeline dependency
bionext_ver="0.7.1"
tool_desc_ver="1.4"

# Check Nextflow
if (command -v nextflow > /dev/null 2>&1); then {
	echo "Nextflow installed"
} else {
  		echo "Please install Nextflow by using:
curl -s https://get.nextflow.io | bash 
"
}
fi
# Check BioNextflow
if [ -d "lib" ]; then {
	echo "BioNextflow installed"
} else {
  		echo "BioNextflow not installed. Installing it..."
		curl -k -L https://github.com/biocorecrg/BioNextflow/archive/${bionext_ver}.tar.gz > ${bionext_ver}.tar.gz
		tar -zvxf ${bionext_ver}.tar.gz
		mv "BioNextflow-${bionext_ver}/lib" .
		rm ${bionext_ver}.tar.gz
		rm -fr BioNextflow-${bionext_ver}
}
fi
#Check Tool description for multiQC
if [ -e "conf_tools.txt" ]; then {
        echo "make_tool_desc_for_multiqc installed"
} else {
                echo "make_tool_desc_for_multiqc not installed. Installing it..."
                curl -k -L https://github.com/CRG-CNAG/make_tool_desc_for_multiqc/archive/${tool_desc_ver}.tar.gz >  ${tool_desc_ver}.tar.gz
                tar -zvxf ${tool_desc_ver}.tar.gz
		mv make_tool_desc_for_multiqc-${tool_desc_ver}/conf_tools.txt .
                rm -fr make_tool_desc_for_multiqc-${tool_desc_ver} 
		rm -f ${tool_desc_ver}.tar.gz
}
fi
