function populateFamilyAndHaploMap(text_unformatted, flow = false){
    //console.log("::"+(flow?"Flow":"Chr")+" --- start");

    var lines = text_unformatted.split('\n');

    var tmp = {
	_fam : null,
	_perc_array : [], // horizontal percs
	_alleles_array : [] // vertical haplos [ [[],[]] , [[],[]] ]
    }



    function flushTmpData(tmp){
	// Finish populating alleles and insert percs
	if (tmp._perc_array.length > 0)
	{
	    if (tmp._perc_array.length !== tmp._alleles_array.length){
		console.log("Length mismatch");
		throw new Error("");
	    }


	    for (var tpa=0; tpa < tmp._perc_array.length; tpa ++)
	    {
		var perc_alleles = tmp._alleles_array[tpa],
		    perc = tmp._perc_array[tpa];


		if (flow){ // flow relies on prior perc existence
		    var perc = familyMapOps.getPerc(perc.id, tmp._fam);
		    perc.insertFlowData( perc_alleles[0] );
		    perc.insertFlowData( perc_alleles[1] );
		    //console.log("INSERTING FLOW", perc.id, perc.haplo_data[0].flow);
		}
		else {
		    perc.insertHaploData( perc_alleles[0] )
		    perc.insertHaploData( perc_alleles[1] )
		    familyMapOps.insertPerc(perc, Number(tmp._fam));
		}
	    }

	    tmp._perc_array = [];
	    tmp._alleles_array = [];
	}
    }




    // Populate Marker
    for (var l=0; l < lines.length; l++)
    {
	var line = lines[l];

	if (line.startsWith("FAMILY")){

	    flushTmpData(tmp);

	    var fid = line.split(/\s+/)[1]
	    tmp._fam = Number(fid);
	    continue
	}


	if (tmp._perc_array.length === 0){
	    // Hunt for names
	    if (line.indexOf('(')!==-1 && line.indexOf(')')!==-1)
	    {
		var people = line.trim().split(/\s{2,}/);

		for (var p=0; p < people.length; p++)
		{
		    var perc = people[p].split(" ");

		    var id = Number(perc[0]),
			parents = perc[1].split("(")[1].split(")")[0];

		    var mother_id = 0,
			father_id = 0;

		    if (parents !== "F"){
			parents = parents.split(",").map(x => Number(x));
			
			mother_id = parents[0];
			father_id = parents[1];
		    }

		    // Gender's and Affecteds are unknown.
		    // Gender's can be inferred, but affectation needs a ped file
		    var newp = new Person(id, 0, 0, mother_id, father_id);
		    tmp._perc_array.push(newp);
		    tmp._alleles_array.push([ [],[] ])
		}
	    }
	    continue;
	}

	var trimmed = line.trim();
	if (trimmed.length == 0){
	    flushTmpData(tmp);
	    continue;
	}

	//Allele lines
	var multiple_alleles = trimmed.split(/\s{3,}/);

	if (multiple_alleles.length  !== tmp._perc_array.length){
	    console.log(trimmed, multiple_alleles, tmp._perc_array)
	    throw new Error("Num alleles and num percs do not align");
	}

	for (var a=0; a < multiple_alleles.length; a++)
	{
	    var alleles = multiple_alleles[a],
		left_right = null;

	    if (!flow){
		// We ignore all types of phasing and for 
		// ambiguously marker alleles "A", we pick the
		// first (this holds consistent for inherited).
		//var left_right = alleles.split(/\s[+:|\\/]\s/)

		left_right = alleles.split(/\s[^\d]\s/)
		    .map(x=> Number(
			x.split(",")[0]
			    .replace("A","")
			//.replace("?","9")
		    ));
	    }
	    else {
		left_right = alleles.split(/\s[^\d]\s/);
		FlowResolver.convertGroupsToFamilySpecific(left_right, tmp._fam);
	    }

	    tmp._alleles_array[a][0].push(left_right[0]);
	    tmp._alleles_array[a][1].push(left_right[1]);
	}
    }
    flushTmpData(tmp);
}
