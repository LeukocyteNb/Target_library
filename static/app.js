const targetGrid = document.getElementById('targetGrid');
const loadLatestBtn = document.getElementById('loadLatestBtn');
const targetStatusEl = document.getElementById('targetStatus');

const designStatusEl = document.getElementById('designStatus');
const designResultEl = document.getElementById('designResult');
const designMetaEl = document.getElementById('designMeta');
const binderSeqEl = document.getElementById('binderSeq');
const bindingNoteEl = document.getElementById('bindingNote');
const designFilesEl = document.getElementById('designFiles');
const viewerEl = document.getElementById('viewer3d');

let viewer = null;
let selectedTarget = '';

function setTargetStatus(message, error = false) {
  targetStatusEl.textContent = message;
  targetStatusEl.style.color = error ? '#a23333' : '#5f6f63';
}

function setDesignStatus(message, error = false) {
  designStatusEl.textContent = message;
  designStatusEl.style.color = error ? '#a23333' : '#5f6f63';
}

function highlightSelectedTarget(receptor) {
  const cards = targetGrid.querySelectorAll('.target-card');
  cards.forEach((card) => {
    card.classList.toggle('active', card.dataset.receptor === receptor);
  });
}

function selectTarget(receptor) {
  selectedTarget = receptor;
  highlightSelectedTarget(receptor);
  designResultEl.classList.add('hidden');
  setDesignStatus(`Selected receptor ${receptor}. You can load latest mini-binder.`);
}

function renderTargetCards(targets) {
  targetGrid.innerHTML = '';

  if (!targets || targets.length === 0) {
    setTargetStatus('No target metadata found. Add entries in data/receptor_structures.json.', true);
    return;
  }

  targets.forEach((target) => {
    const card = document.createElement('button');
    card.type = 'button';
    card.className = 'target-card';
    card.dataset.receptor = target.receptor;
    const tagsHtml = (target.process_tags || [])
      .map((tag) => `<span class="tag-chip">${tag}</span>`)
      .join('');
    card.innerHTML = `
      <span class="target-name">${target.label || target.receptor}</span>
      <span class="target-meta">${target.receptor} | ${target.target_family || 'Unknown'}</span>
      <span class="target-meta">${target.has_design ? 'Latest design available' : 'No local design yet'}</span>
      <span class="target-tags">${tagsHtml}</span>
    `;

    card.addEventListener('click', () => {
      selectTarget(target.receptor);
      setTargetStatus(`Selected ${target.receptor}.`);
    });

    targetGrid.appendChild(card);
  });
}

async function fetchTargetCatalog() {
  setTargetStatus('Loading target catalog...');
  try {
    const response = await fetch('/api/targets');
    const payload = await response.json();
    if (!response.ok) {
      throw new Error(payload.error || 'Failed to load target catalog');
    }
    renderTargetCards(payload.targets || []);
    setTargetStatus(`Loaded ${payload.targets?.length || 0} targets.`);
  } catch (err) {
    setTargetStatus(`Failed to load target catalog: ${err.message}`, true);
  }
}

function ensureViewer() {
  if (!window.$3Dmol) {
    setDesignStatus('3Dmol library not loaded. Please refresh or check access to 3Dmol.org.', true);
    return null;
  }
  if (!viewer) {
    viewer = $3Dmol.createViewer(viewerEl, { backgroundColor: '#f6faf5' });
  }
  return viewer;
}

function render3D(complexPdb) {
  const v = ensureViewer();
  if (!v || !complexPdb) {
    setDesignStatus('3D viewer is unavailable in this browser session.', true);
    return;
  }

  try {
    v.clear();
    v.addModel(complexPdb, 'pdb');
    v.setStyle({ chain: 'A' }, { cartoon: { color: '#6f9b73' }, line: { color: '#6f9b73', linewidth: 2.0 } });
    v.setStyle({ chain: 'B' }, { stick: { color: '#c96f3a', radius: 0.28 }, line: { color: '#c96f3a', linewidth: 2.2 } });
    v.zoomTo();
    v.resize();
    v.render();
  } catch (err) {
    setDesignStatus(`3D render failed: ${err.message}`, true);
  }
}

function renderDesignResult(payload) {
  designResultEl.classList.remove('hidden');

  const binding = payload.binding || {};
  const hotspots = (binding.target_hotspot_residues || []).join(', ') || 'n/a';

  designMetaEl.textContent = `Design ID: ${payload.design_id} | Receptor: ${payload.receptor} | Engine: ${payload.engine}`;
  binderSeqEl.textContent = payload.binder_sequence || 'N/A';
  bindingNoteEl.textContent = `${binding.binding_interface_note || 'No binding note.'} Hotspots: ${hotspots}. Estimated distance: ${binding.distance_estimate_angstrom || 'n/a'} A.`;

  const files = payload.files || {};
  designFilesEl.innerHTML = `
    <a href="${files.complex_pdb}" target="_blank" rel="noopener noreferrer">Download complex PDB</a>
    <a href="${files.binder_pdb}" target="_blank" rel="noopener noreferrer">Download binder PDB</a>
    <a href="${files.report_json}" target="_blank" rel="noopener noreferrer">Download report JSON</a>
  `;

  render3D(payload.complex_pdb_text || '');
}

async function loadLatestForTarget() {
  if (!selectedTarget) {
    setTargetStatus('Select a target first.', true);
    return;
  }

  setTargetStatus(`Loading latest design for ${selectedTarget}...`);
  setDesignStatus(`Loading latest design for ${selectedTarget}...`);

  try {
    const params = new URLSearchParams({ receptor: selectedTarget });
    const response = await fetch(`/api/target_design?${params.toString()}`);
    const payload = await response.json();

    if (!response.ok) {
      if (response.status === 404) {
        throw new Error(`No local design found for ${selectedTarget}.`);
      }
      throw new Error(payload.error || 'Failed to load latest design');
    }

    renderDesignResult(payload);
    setDesignStatus(`Loaded latest design: ${payload.design_id}`);
    setTargetStatus(`Loaded latest design for ${selectedTarget}.`);
  } catch (err) {
    setTargetStatus(err.message, true);
    setDesignStatus(err.message, true);
  }
}

loadLatestBtn.addEventListener('click', loadLatestForTarget);
fetchTargetCatalog();
