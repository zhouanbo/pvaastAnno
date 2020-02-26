const parse = require('csv-parse/lib/sync')
const stringify = require('csv-stringify/lib/sync')
const fs = require('fs')
const execSync = require('child_process').execSync;

const firebase = require("firebase")
require("firebase/firestore")
firebase.initializeApp({
    apiKey: 'AIzaSyCgLHGZVLRw_V4pclj7voTwp2a4GbnXPa0',
    authDomain: 'geneanno.firebaseapp.com',
    projectId: 'geneanno'
});

const adhdgene = require("../res/adhdgene").adhdgene
const sfari = require("../res/sfari").sfari
const hom_convert = require("../res/hom_conversion").hom_convert
const ncbi_convert = require("../res/ncbi_conversion").ncbi_convert
const refseq_convert = JSON.parse(require("../res/refseq_conversion").refseq_convert)
const synonym_convert = require("../res/synonym_conversion").synonym_convert


const main = async () => {
    
    //read csv file
    console.log(`Reading csv...`)
    let genes = []
    let file = fs.readFileSync(process.argv[2], 'utf-8')
    const records = parse(file, {relax_column_count: true})
    records.map((line,i) => {
        if(line[1] && i>1) genes.push(line[1])
    })

    //search for synonym
    genes = genes.map(gene => {
        if (synonym_convert[gene]) {
        return synonym_convert[gene]
        } else {
        return gene
        }
    })

    //process genes
    let table = {}, snapshot, r
    for(let gene of genes) {
        console.log(`Processing...(${gene})`)
        table[gene] = {}

        //Refseq
        snapshot = await firebase.firestore().collection('refseq').where('name', '==', gene.toUpperCase()).get()
        r = []

        snapshot.forEach(doc => {
        let data = doc.data()
        let chr = refseq_convert.find(e => e.refseq_id == data.contig)
        if (!chr) return
        r.push({
            name: data.name,
            description: decodeURIComponent(data.description),
            chr: chr.ucsc_id,
            contig: data.contig,
            start: data.start,
            end: data.end,
            type: data.biotype,
            strand: data.strand,
            synonym: data.synonym
        })
        })
        if (r.length > 0)
        table[gene] = Object.assign(table[gene], r.sort((a,b) => a.contig > b.contig ? 1 : -1)[0])


        //ExAC
        snapshot = await firebase.firestore().collection('exac').where('gene', '==', gene.toUpperCase()).orderBy('pLI', 'desc').limit(1).get()
        r = []
        snapshot.forEach(doc => {
        let data = doc.data()
        table[gene] = Object.assign(table[gene], {
            pLI: data.pLI
        })
        })


        //SFARI
        table[gene] = Object.assign(table[gene], {
        sfari: sfari[gene.toUpperCase()] ? sfari[gene.toUpperCase()] : 0
        })


        //ADHDgene
        table[gene] = Object.assign(table[gene], {
        adhdgene: adhdgene[gene.toUpperCase()] ? 1 : 0
        })


        //DISEASE
        snapshot = await firebase.firestore().collection('disease_textmining').where('SYMBOL', '==', gene.toUpperCase()).orderBy('Confidence', 'desc').get()
        r = []
        snapshot.forEach(doc => {
        let data = doc.data()
        r.push({
            name: decodeURIComponent(data.NAME),
            value: data.Confidence,
        })
        })
        if (r.length > 0)
        table[gene] = Object.assign(table[gene], {tm: r})

        snapshot = await firebase.firestore().collection('disease_experiments').where('SYMBOL', '==', gene.toUpperCase()).orderBy('Confidence', 'desc').get()
        r = []
        snapshot.forEach(doc => {
        let data = doc.data()
        r.push({
            name: decodeURIComponent(data.NAME),
            value: data.Confidence,
        })
        })
        if (r.length > 0)
        table[gene] = Object.assign(table[gene], {exp: r})

        snapshot = await firebase.firestore().collection('disease_knowledge').where('SYMBOL', '==', gene.toUpperCase()).orderBy('Confidence', 'desc').get()
        r = []
        snapshot.forEach(doc => {
        let data = doc.data()
        r.push({
            name: decodeURIComponent(data.NAME),
            value: data.Confidence,
        })
        })
        if (r.length > 0)
        table[gene] = Object.assign(table[gene], {kl: r})


        //IMPC
        if (hom_convert[gene.toUpperCase()]){
        snapshot = await firebase.firestore().collection('impc').where('marker_symbol', '==', hom_convert[gene.toUpperCase()]).get()
        r = []
        snapshot.forEach(doc => {
            let data = doc.data()
            if (data.mp_term_name)
            r.push({
                name: decodeURIComponent(data.mp_term_name),
                value: data.p_value == 0 ? 0.0000000001 : data.p_value,
            })
        })
        r = r.filter((e, idx) =>
            idx === r.findIndex((t) => (
            t.name === e.name
            ))
        ).map(e => {
            e.value = -Math.log10(e.value)
            e.tvalue = e.value
            if(e.value>6) e.value = 6
            return e
        })
        if (r.length > 0)
            table[gene] = Object.assign(table[gene], {ho: r.sort((a,b) => a.tvalue < b.tvalue ? 1 : -1)})
        }


        //HumanBase
        if (ncbi_convert[gene.toUpperCase()]) {
        await fetch(`https://us-central1-geneanno.cloudfunctions.net/humanbase/?id=${ncbi_convert[gene.toUpperCase()]}`).then((response) => {
            return response.json();
        })
        .then((myJson) => {
            let data = myJson.sort((a,b) => a.score < b.score ? 1 : -1)
            r = data.map(e => {
            e.name = decodeURIComponent(e.term.title)
            e.value = e.score
            return e
            })
            if (r.length > 0)
            table[gene] = Object.assign(table[gene], {hb: r});
        });
        }


        //GTEx
        snapshot = await firebase.firestore().collection('gtex').where('Description', '==', gene.toUpperCase()).limit(1).get()
        r = []
        snapshot.forEach(doc => {
        let data = doc.data()
        delete data.Name
        delete data.Description
        Object.keys(data).map(e => {
            r.push({
            name: decodeURIComponent(e),
            value: data[e]
            })
        })
        })
        if (r.length > 0)
        table[gene] = Object.assign(table[gene], {gt: r.sort((a,b) => a.value < b.value ? 1 : -1)})


        //Allen
        snapshot = await firebase.firestore().collection('allen').where('name', '==', gene.toUpperCase()).limit(1).get()
        r = []
        snapshot.forEach(doc => {
        let data = doc.data()
        delete data.name
        Object.keys(data).map(e => {
            r.push({
            name: decodeURIComponent(e),
            value: data[e]
            })
        })
        })
        if (r.length > 0)
        table[gene] = Object.assign(table[gene], {al: r.sort((a,b) => a.value < b.value ? 1 : -1)})
        
    }//end of gene for loop

  
    //variants AF
    const findgene = (n) => {
        while(!records[n][1]){n--}
        return n
    }
    records.map((line, i) => {
        if (i>1) {        
            let chr = records[findgene(i)][5].replace('chr', '')
            let pos = line[15]
            let geno = line[18].split("|")
            let ref = line[17].split("|")[0]
            let alt = ''
            for (e of geno) {
                if (e.includes(':') && !e.includes('^')) {
                    let arr = e.split(":")
                    if (arr[0] !== ref) alt = arr[0]
                    if (arr[1] !== ref) alt = arr[1]
                    break
                } 
            }
            fs.appendFileSync('tmpVar.txt', chr+'\t'+pos+'\t'+(Number(pos)+ref.length-1)+'\t'+ref+'\t'+alt+'\n')
        }
    })
    console.log("Searching gnomad...")
    let gnomad_grep = execSync(`grep -f tmpVar.txt ../res/hg19_gnomad_exome.txt`, { encoding: 'utf-8' }).split('\n');
    console.log("Searching exac_nonpsych...")
    let exac_grep = execSync(`grep -f tmpVar.txt ../res/hg19_exac03nonpsych.txt`, { encoding: 'utf-8' }).split('\n');
    console.log("Searching PolyPhen...")
    let pp_grep = execSync(`grep -f tmpVar.txt ../res/hg19_ljb2_pp2hdiv.txt`, { encoding: 'utf-8' }).split('\n');
    console.log("Searching SIFT...")
    let sift_grep = execSync(`grep -f tmpVar.txt ../res/hg19_ljb23_sift.txt`, { encoding: 'utf-8' }).split('\n');
    gnomad_grep = gnomad_grep.map(l => l.split('\t'))
    exac_grep = exac_grep.map(l => l.split('\t'))
    pp_grep = pp_grep.map(l => l.split('\t'))
    sift_grep = sift_grep.map(l => l.split('\t'))
    fs.unlinkSync('tmpVar.txt')
    
    
    //generate csv
    col = ["StandardName","Start","End","Description","Type","Synonym","pLI","ASD","ADHD","Textmining","Experiments","Knowledge","MousePheno","TissueML","GTEx","BrainSpan"]
    col2 = col.map(() => "")
    col3 = ["gnomad_AF", "ExAC_nonpsych_AF", "PolyPhen", "SIFT"]

    records[1].splice(11, 0, ...col)
    records[1].splice(19+col.length, 0, ...col3)
    records.map((line, i) => {
        if(line[1] && i>1) {
            let gene = synonym_convert[line[1]] || line[1]
            let csv_array = []
            csv_array.push(
                table[gene].name, 
                table[gene].start || 'NA', 
                table[gene].end || 'NA', 
                table[gene].description || 'NA',
                table[gene].type || 'NA', 
                table[gene].synonym?table[gene].synonym:'NA', 
                table[gene].pLI?table[gene].pLI.toFixed(2):'NA', 
                table[gene].sfari || '0', 
                table[gene].adhdgene || '0', 
                table[gene].tm?table[gene].tm.map(e=>(e.tvalue?e.tvalue:e.value).toFixed(2)+":"+e.name).join('|'):'NA', 
                table[gene].exp?table[gene].exp.map(e=>(e.tvalue?e.tvalue:e.value).toFixed(2)+":"+e.name).join('|'):'NA', 
                table[gene].kl?table[gene].kl.map(e=>(e.tvalue?e.tvalue:e.value).toFixed(2)+":"+e.name).join('|'):'NA',
                table[gene].ho?table[gene].ho.map(e=>(e.tvalue?e.tvalue:e.value).toFixed(2)+":"+e.name).join('|'):'NA',
                table[gene].hb?table[gene].hb.map(e=>(e.tvalue?e.tvalue:e.value).toFixed(2)+":"+e.name).join('|'):'NA',
                table[gene].gt?table[gene].gt.map(e=>(e.tvalue?e.tvalue:e.value).toFixed(2)+":"+e.name).join('|'):'NA',
                table[gene].al?table[gene].al.map(e=>(e.tvalue?e.tvalue:e.value).toFixed(2)+":"+e.name).join('|'):'NA'
            )
            records[i].splice(11, 0, ...csv_array)
        } else if (!line[1]) {
            records[i].splice(11, 0, ...col2)
        }
        if(i>1){
            let csv_array = []
            let chr = records[findgene(i)][5].replace('chr', '')
            let pos = line[15+col.length]
            let geno = line[18+col.length].split("|")
            let ref = ''
            let alt = ''
            for (e of geno) {
                if (e.includes(':') && !e.includes('^')) {
                    let arr = e.split(":")
                    ref = arr[0]
                    alt = arr[1]
                    break
                } 
            }
            let gnomad = gnomad_grep.filter(line => line[0]==chr && line[1]==pos && line[3]==ref && line[4]==alt)[0]
            let exac = exac_grep.filter(line => line[0]==chr && line[1]==pos && line[3]==ref && line[4]==alt)[0]
            let pp = pp_grep.filter(line => line[0]==chr && line[1]==pos && line[3]==ref && line[4]==alt)[0]
            let sift = sift_grep.filter(line => line[0]==chr && line[1]==pos && line[3]==ref && line[4]==alt)[0]
            gnomad = gnomad?gnomad[5]:'NA'
            exac = exac?exac[5]:'NA'
            pp = pp?pp[6]:'NA'
            sift = sift?sift[7]:'NA'
            csv_array.push(
                gnomad,
                exac,
                pp,
                sift
            )
            records[i].splice(19+col.length, 0, ...csv_array)
        }
    })

    
    fs.writeFileSync("output.csv", stringify(records))

}

main().then(() => {
    process.exit()
}).catch(e => {
    console.error(e)
})
