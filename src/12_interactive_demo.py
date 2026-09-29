"""Day 6: build a self-contained interactive path explorer (docs/interactive/index.html).

Reads the tables written by 09/10 (candidates_<slug>_triaged.csv, explain_summary_<slug>.csv, paths_<slug>.csv)
and embeds them in one static HTML page (no server, no CDN): pick a disease and a candidate, see the supporting
paths as a layered graph, toggle metapaths, change how many paths are drawn, hover/click a node to see the paths
through it, and read the gene-ablation and triage notes. Candidates excluded by triage are listed (greyed) with
the reason. Paths through disease-named pathways are hidden (see 07b).

Usage:
    python3 src/12_interactive_demo.py
Then open docs/interactive/index.html (or publish the folder with GitHub Pages).
"""
from __future__ import annotations

import json
import re

import numpy as np
import pandas as pd

from dwpc import build_matrices
from utils import OUT_TAB, ROOT, load_processed

PARSE = re.compile(r" -binds-> | <-assoc- | -interacts- | -in-> | <-in- ")
DISEASES = {"alzheimer": "Alzheimer's disease", "parkinson": "Parkinson's disease"}
OUT = ROOT / "docs" / "interactive"
EXCLUDED = {"endogenous", "diagnostic", "cytotoxic_or_metal", "risk_factor", "wrong_direction_likely"}


def num(x, nd=3):
    return None if pd.isna(x) else round(float(x), nd)


def build_disease(slug: str, hide: set) -> dict | None:
    cpath = OUT_TAB / f"candidates_{slug}_triaged.csv"
    spath = OUT_TAB / f"explain_summary_{slug}.csv"
    ppath = OUT_TAB / f"paths_{slug}.csv"
    if not (cpath.exists() and spath.exists() and ppath.exists()):
        print(f"  skipping {slug}: missing one of {cpath.name}, {spath.name}, {ppath.name}")
        return None
    cand = pd.read_csv(cpath).fillna({"category": "", "note": ""})
    summ = pd.read_csv(spath).set_index("compound")
    paths = pd.read_csv(ppath)
    out = []
    for rank, r in enumerate(cand.itertuples(index=False), start=1):
        d = {"compound": r.compound, "order": rank, "consensus": num(r.consensus),
             "pct_dwpc": num(r.pct_dwpc), "pct_z": num(r.pct_z), "n_bound_genes": int(r.n_bound_genes),
             "category": r.category, "note": r.note, "shortlisted": bool(r.shortlisted),
             "excluded": r.category in EXCLUDED, "paths": []}
        if r.shortlisted and r.compound in summ.index:
            s = summ.loc[r.compound]
            d.update({
                "tier": "multi-gene" if s["rank_no_top1"] <= 0.05 * s["n_ranked"] else "single-gene-dependent",
                "top_gene": s["top_gene"], "top_gene_share": num(s["top_gene_share"], 2),
                "top3_genes": s["top3_genes"], "dwpc_rank": int(s["dwpc_rank"]),
                "rank_no_top1": int(s["rank_no_top1"]), "rank_no_top3": int(s["rank_no_top3"]),
                "n_ranked": int(s["n_ranked"]),
                "n_paths": {"CbGaD": int(s["n_paths_CbGaD"]), "CbGiGaD": int(s["n_paths_CbGiGaD"]),
                            "CbGpPWpGaD": int(s["n_paths_CbGpPWpGaD"])}})
            sub = paths[paths["compound"] == r.compound]
            for p in sub.itertuples(index=False):
                parts = PARSE.split(p.path)
                if p.metapath == "CbGpPWpGaD" and parts[2] in hide:
                    continue
                d["paths"].append({"mp": p.metapath, "share": num(p.share_in_metapath, 4),
                                   "w": float(f"{p.weight:.4g}"), "n": parts})
        out.append(d)
    return {"name": DISEASES[slug], "candidates": out}


def main() -> None:
    nodes, edges = load_processed()
    hide = set(build_matrices(nodes, edges, drop_disease_pathways=True).dropped_pathways or [])
    data = {}
    for slug in DISEASES:
        d = build_disease(slug, hide)
        if d:
            data[slug] = d
    if not data:
        raise SystemExit("No inputs found - run 09 and 10 first.")
    OUT.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(data, separators=(",", ":")).replace("</", "<\\/")
    (OUT / "index.html").write_text(HTML.replace("__DATA__", payload), encoding="utf-8")
    n = sum(len(v["candidates"]) for v in data.values())
    print(f"Wrote {OUT / 'index.html'}  ({len(data)} diseases, {n} candidates)")


HTML = r"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Explainable drug repurposing paths - Hetionet (AD / PD)</title>
<style>
:root{--bg:#fafafa;--fg:#222;--muted:#666;--card:#fff;--line:#ddd;--acc:#1f77b4}
@media (prefers-color-scheme: dark){:root{--bg:#161616;--fg:#e8e8e8;--muted:#9a9a9a;--card:#1f1f1f;--line:#333}}
*{box-sizing:border-box}body{margin:0;font:14px/1.45 system-ui,-apple-system,Segoe UI,Roboto,sans-serif;background:var(--bg);color:var(--fg)}
header{padding:14px 18px;border-bottom:1px solid var(--line)}h1{font-size:18px;margin:0 0 4px}header p{margin:0;color:var(--muted);font-size:13px}
main{display:grid;grid-template-columns:290px 1fr 330px;gap:12px;padding:12px;align-items:start}
@media (max-width:1150px){main{grid-template-columns:1fr}}
.card{background:var(--card);border:1px solid var(--line);border-radius:8px;padding:10px}
#list{max-height:78vh;overflow:auto}
.item{padding:6px 8px;border-radius:6px;cursor:pointer;display:flex;justify-content:space-between;gap:6px;align-items:center}
.item:hover{background:rgba(31,119,180,.12)}.item.sel{background:rgba(31,119,180,.22)}.item.dim{opacity:.55}
.badge{font-size:10.5px;padding:1px 6px;border-radius:9px;border:1px solid var(--line);white-space:nowrap;color:var(--muted)}
.b-ok{border-color:#2ca02c;color:#2ca02c}.b-bad{border-color:#d62728;color:#d62728}.b-ex{border-color:#999}
select,input[type=range]{width:100%}label{display:block;margin:4px 0;color:var(--muted)}
svg{width:100%;height:auto;background:var(--card)}
.node{cursor:pointer}.node text{font-size:11px;fill:var(--fg);pointer-events:none}
.edge{fill:none;stroke-linecap:round}.dimmed{opacity:.12}
#tip{position:fixed;pointer-events:none;background:#000c;color:#fff;padding:5px 8px;border-radius:5px;font-size:12px;display:none;max-width:320px;z-index:5}
.kv{display:grid;grid-template-columns:auto 1fr;gap:2px 10px;margin:6px 0}.kv b{color:var(--muted);font-weight:500}
ol{padding-left:18px;margin:6px 0}li{margin:3px 0;font-size:12.5px}.legend span{display:inline-flex;align-items:center;margin-right:10px;font-size:12px;color:var(--muted)}
.legend i{width:10px;height:10px;border-radius:50%;display:inline-block;margin-right:4px;border:1px solid #0006}
.note{font-size:12px;color:var(--muted);margin-top:8px}
</style></head><body>
<header><h1>Explainable drug-repurposing paths (Hetionet v1.0)</h1>
<p>Precomputed graph evidence for the top consensus candidates. Paths show why the graph suggests a compound - not mechanism, and not direction of effect (binding does not distinguish agonist from antagonist).</p></header>
<main>
<section class="card"><label>Disease<select id="dis"></select></label>
<div id="list"></div></section>
<section class="card"><div id="head"></div>
<div class="legend"><span><i style="background:#4c78a8"></i>compound</span><span><i style="background:#9ecae9"></i>bound gene</span><span><i style="background:#a1d99b"></i>pathway</span><span><i style="background:#fdae6b"></i>disease gene</span><span><i style="background:#e6550d"></i>disease</span></div>
<svg id="svg" viewBox="0 0 1000 560" preserveAspectRatio="xMidYMid meet"></svg>
<div class="note" id="foot"></div></section>
<section class="card"><div id="ctl"></div><div id="side"></div></section>
</main><div id="tip"></div>
<script id="data" type="application/json">__DATA__</script>
<script>
const DATA = JSON.parse(document.getElementById('data').textContent);
const MP = {
 CbGaD:{roles:['compound','dgene','disease'],rels:['binds','assoc'],label:'CbGaD - binds a disease gene'},
 CbGiGaD:{roles:['compound','bgene','dgene','disease'],rels:['binds','interacts','assoc'],label:'CbGiGaD - binds a gene that interacts with a disease gene'},
 CbGpPWpGaD:{roles:['compound','bgene','pathway','dgene','disease'],rels:['binds','in','in','assoc'],label:'CbGpPWpGaD - binds a gene sharing a pathway with a disease gene'}
};
const LAYER={compound:0,bgene:1,pathway:2,dgene:3,disease:4};
const NCOL={compound:'#4c78a8',bgene:'#9ecae9',pathway:'#a1d99b',dgene:'#fdae6b',disease:'#e6550d'};
const ECOL={binds:'#1f77b4',interacts:'#8a8a8a',in:'#2ca02c',assoc:'#d62728'};
const state={dis:Object.keys(DATA)[0],cand:null,perMp:3,on:{CbGaD:true,CbGiGaD:true,CbGpPWpGaD:true},pin:null};
const $=id=>document.getElementById(id), NS='http://www.w3.org/2000/svg';
const tip=$('tip');
function el(tag,attrs,parent){const e=document.createElementNS(NS,tag);for(const k in attrs)e.setAttribute(k,attrs[k]);if(parent)parent.appendChild(e);return e}
function esc(s){return String(s).replace(/[&<>"]/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;'}[c]))}
function fmt(x,n=2){return x==null?'-':Number(x).toFixed(n)}

function init(){
  const sel=$('dis');
  for(const k in DATA){const o=document.createElement('option');o.value=k;o.textContent=DATA[k].name;sel.appendChild(o)}
  sel.onchange=()=>{state.dis=sel.value;state.cand=firstShort();state.pin=null;renderAll()};
  state.cand=firstShort();renderAll();
}
function cands(){return DATA[state.dis].candidates}
function firstShort(){const c=cands().find(x=>x.shortlisted);return c?c.compound:cands()[0].compound}
function cur(){return cands().find(x=>x.compound===state.cand)}

function renderAll(){renderList();renderCtl();renderCand()}
function renderList(){
  const box=$('list');box.innerHTML='';
  cands().forEach(c=>{
    const d=document.createElement('div');d.className='item'+(c.compound===state.cand?' sel':'')+(c.shortlisted?'':' dim');
    let b='';
    if(c.shortlisted)b=`<span class="badge ${c.tier==='multi-gene'?'b-ok':'b-bad'}">${c.tier==='multi-gene'?'multi-gene':'1-gene'}</span>`;
    else if(c.category)b=`<span class="badge ${c.excluded?'b-ex':''}">${esc(c.category.replace(/_/g,' '))}</span>`;
    else b=`<span class="badge">not reviewed</span>`;
    d.innerHTML=`<span>${c.order}. ${esc(c.compound)}</span>${b}`;
    d.onclick=()=>{state.cand=c.compound;state.pin=null;renderList();renderCand()};
    box.appendChild(d);
  });
}
function renderCtl(){
  const c=$('ctl');
  c.innerHTML=`<b>Display</b>`+Object.keys(MP).map(k=>`<label><input type="checkbox" data-mp="${k}" ${state.on[k]?'checked':''}> ${MP[k].label}</label>`).join('')+
  `<label>Paths per metapath: <b id="pm">${state.perMp}</b><input id="per" type="range" min="1" max="10" value="${state.perMp}"></label>`;
  c.querySelectorAll('input[type=checkbox]').forEach(i=>i.onchange=()=>{state.on[i.dataset.mp]=i.checked;state.pin=null;renderCand()});
  $('per').oninput=e=>{state.perMp=+e.target.value;$('pm').textContent=state.perMp;state.pin=null;renderCand()};
}

function selectedPaths(c){
  const out=[];
  for(const mp of Object.keys(MP)){if(!state.on[mp])continue;
    out.push(...c.paths.filter(p=>p.mp===mp).sort((a,b)=>b.share-a.share).slice(0,state.perMp))}
  return out;
}
function buildGraph(paths){
  const nodes=new Map(),edges=new Map();
  paths.forEach((p,pi)=>{const m=MP[p.mp];
    const keys=p.n.map((lab,i)=>lab+'|'+m.roles[i]);
    keys.forEach((k,i)=>{if(!nodes.has(k))nodes.set(k,{key:k,label:p.n[i],role:m.roles[i],layer:LAYER[m.roles[i]]})});
    for(let i=0;i<keys.length-1;i++){const ek=keys[i]+'>'+keys[i+1];
      if(!edges.has(ek))edges.set(ek,{a:keys[i],b:keys[i+1],rel:m.rels[i],share:0,paths:[]});
      const e=edges.get(ek);e.share+=p.share;e.paths.push(pi)}
  });
  return{nodes,edges};
}
function layout(g,W,H){
  const layers=[[],[],[],[],[]];g.nodes.forEach(n=>layers[n.layer].push(n));
  const y=new Map();
  const assign=()=>layers.forEach(L=>{L.forEach((n,i)=>y.set(n.key,i))});
  const nbr=(n,dir)=>[...g.edges.values()].filter(e=>dir<0?e.b===n.key:e.a===n.key).map(e=>dir<0?e.a:e.b);
  assign();
  const bary=(n,dir)=>{const ns=nbr(n,dir);if(!ns.length)return y.get(n.key);
    return ns.reduce((s,k)=>s+ (y.get(k)+ (maxN-layers[g.nodes.get(k).layer].length)/2),0)/ns.length - (maxN-layers[n.layer].length)/2};
  let maxN=Math.max(...layers.map(L=>L.length),1);
  for(let it=0;it<4;it++){
    for(let l=1;l<5;l++){const b=new Map(layers[l].map(n=>[n.key,bary(n,-1)]));layers[l].sort((p,q)=>b.get(p.key)-b.get(q.key));assign()}
    for(let l=3;l>=0;l--){const b=new Map(layers[l].map(n=>[n.key,bary(n,1)]));layers[l].sort((p,q)=>b.get(p.key)-b.get(q.key));assign()}
  }
  const padX=70,padY=36,sp=Math.min(46,(H-2*padY)/Math.max(maxN-1,1));
  layers.forEach((L,l)=>L.forEach((n,i)=>{n.x=padX+l*(W-2*padX)/4;n.y=H/2+(i-(L.length-1)/2)*sp}));
  return maxN;
}
function renderCand(){
  const c=cur(),svg=$('svg');svg.innerHTML='';
  $('head').innerHTML=`<b style="font-size:16px">${esc(c.compound)}</b> <span class="badge">${esc(DATA[state.dis].name)}</span> `+
    (c.shortlisted?`<span class="badge ${c.tier==='multi-gene'?'b-ok':'b-bad'}">${c.tier}</span>`:`<span class="badge b-ex">${esc((c.category||'not explained').replace(/_/g,' '))}</span>`);
  if(!c.shortlisted){
    $('foot').textContent='';
    $('side').innerHTML=sideHTML(c,[])+`<p class="note">${c.excluded?`Excluded from the explained shortlist by triage (${esc(c.category.replace(/_/g,' '))}).`:'Outside the explained shortlist (paths were only computed for the top candidates that pass triage).'}</p>`;
    const t=el('text',{x:500,y:280,'text-anchor':'middle','font-size':16,fill:'#888'},svg);t.textContent=c.excluded?'No paths drawn: excluded by triage':'No paths drawn: outside the explained shortlist';return;
  }
  const paths=selectedPaths(c),g=buildGraph(paths);
  if(!paths.length){$('side').innerHTML=sideHTML(c,paths);return}
  const cnt=[0,0,0,0,0];g.nodes.forEach(n=>cnt[n.layer]++);
  const W=1000,H=Math.max(340,Math.max(...cnt)*50+110);svg.setAttribute('viewBox',`0 0 ${W} ${H}`);layout(g,W,H);
  const emax=Math.max(...[...g.edges.values()].map(e=>e.share));
  const eg=el('g',{},svg),ng=el('g',{},svg),E=[];
  g.edges.forEach(e=>{const a=g.nodes.get(e.a),b=g.nodes.get(e.b),span=Math.abs(a.layer-b.layer);
    const dx=(b.x-a.x)*0.5,bow=span>=3?(b.y<a.y?-70:70):0;
    const p=el('path',{class:'edge',d:`M${a.x},${a.y} C${a.x+dx},${a.y+bow} ${b.x-dx},${b.y+bow} ${b.x},${b.y}`,
      stroke:ECOL[e.rel],'stroke-width':1+6*e.share/emax,opacity:.75},eg);
    E.push({e,p});
    p.addEventListener('mousemove',ev=>showTip(ev,`${g.nodes.get(e.a).label} → ${g.nodes.get(e.b).label}<br>relation: ${e.rel} · ${e.paths.length} path(s)`));
    p.addEventListener('mouseleave',()=>tip.style.display='none');
  });
  const N=[];
  g.nodes.forEach(n=>{const gr=el('g',{class:'node',transform:`translate(${n.x},${n.y})`},ng);
    el('circle',{r:n.role==='compound'||n.role==='disease'?18:12,fill:NCOL[n.role],stroke:'#0008','stroke-width':.8},gr);
    const t=el('text',{x:0,y:n.role==='compound'||n.role==='disease'?32:26,'text-anchor':'middle'},gr);
    t.textContent=n.label.length>22?n.label.slice(0,21)+'…':n.label;N.push({n,gr});
    gr.addEventListener('mouseenter',()=>highlight(n.key));gr.addEventListener('mouseleave',()=>highlight(state.pin));
    gr.addEventListener('click',()=>{state.pin=state.pin===n.key?null:n.key;highlight(state.pin);sideRender()});
    gr.addEventListener('mousemove',ev=>showTip(ev,`${esc(n.label)}<br><i>${n.role}</i>`));
    gr.addEventListener('mouseleave',()=>tip.style.display='none');
  });
  function highlight(key){
    E.forEach(({e,p})=>p.classList.toggle('dimmed',key!=null&&e.a!==key&&e.b!==key));
    N.forEach(({n,gr})=>gr.style.opacity=(key==null||n.key===key||[...g.edges.values()].some(e=>(e.a===key&&e.b===n.key)||(e.b===key&&e.a===n.key)))?1:.25);
  }
  function sideRender(){
    let list=paths;
    if(state.pin)list=paths.filter(p=>{const m=MP[p.mp];return p.n.some((lab,i)=>lab+'|'+m.roles[i]===state.pin)});
    $('side').innerHTML=sideHTML(c,list,state.pin&&g.nodes.get(state.pin));
  }
  highlight(state.pin);sideRender();
  $('foot').textContent=`Showing ${paths.length} of ${c.n_paths.CbGaD+c.n_paths.CbGiGaD+c.n_paths.CbGpPWpGaD} supporting paths (top ${state.perMp} per metapath by share of the metapath's weight); edge width ~ share. Click a node to pin it.`;
}
function sideHTML(c,list,pinned){
  let h=`<div class="kv"><b>Consensus</b><span>${fmt(c.consensus,3)} (DWPC pct ${fmt(c.pct_dwpc)}, z pct ${fmt(c.pct_z)})</span><b>Genes bound</b><span>${c.n_bound_genes}</span>`;
  if(c.shortlisted)h+=`<b>Top driving gene</b><span>${esc(c.top_gene)} (${fmt(c.top_gene_share*100,0)}% of weight)</span>
   <b>DWPC-binding rank</b><span>${c.dwpc_rank} of ${c.n_ranked}</span><b>Without top gene</b><span>${c.rank_no_top1}</span><b>Without top 3</b><span>${c.rank_no_top3} <small>(${esc(c.top3_genes)})</small></span>
   <b>Paths (total)</b><span>${c.n_paths.CbGaD} / ${c.n_paths.CbGiGaD} / ${c.n_paths.CbGpPWpGaD}</span>`;
  h+=`</div>`;
  if(c.category||c.note)h+=`<p class="note"><b>Triage${c.category?': '+esc(c.category.replace(/_/g,' ')):''}</b> - ${esc(c.note)}</p>`;
  if(pinned)h+=`<p><b>Paths through ${esc(pinned.label)}</b> (${list.length})</p>`;
  else if(list.length)h+=`<p><b>Displayed paths</b></p>`;
  if(list.length)h+='<ol>'+list.slice().sort((a,b)=>b.share-a.share).map(p=>`<li><small>${p.mp} · ${fmt(p.share*100,0)}%</small><br>${p.n.map(esc).join(' → ')}</li>`).join('')+'</ol>';
  return h;
}
function showTip(ev,html){tip.innerHTML=html;tip.style.display='block';tip.style.left=(ev.clientX+12)+'px';tip.style.top=(ev.clientY+12)+'px'}
init();
</script></body></html>
"""

if __name__ == "__main__":
    main()
