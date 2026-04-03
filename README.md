# Immune receptor ranking + cloud mini-binder workflow

## Start

```bash
python app.py
```

Open: `http://127.0.0.1:8000`

## Deploy to Railway (Option 2)

1. Push this project to a GitHub repository.
2. In Railway, click `New Project` -> `Deploy from GitHub repo`.
3. Select this repository.
4. In service settings, set `Start Command` to:

```bash
python app.py
```

5. Railway will inject `PORT` automatically. The app now listens on `0.0.0.0:$PORT`.
6. After deploy, open the generated Railway domain and share it with friends.

### Optional env vars

- `HOST` (default `0.0.0.0`)
- `PORT` (set by Railway)

## What changed

The app now supports 2 execution modes for mini-binder design:

- `Cloud (Recommended)`
- `Auto / Local`

Use the selector in the design panel.

## Cloud mode setup (simple)

1. Open `data/cloud_backend.example.json`.
2. Copy it to `data/cloud_backend.json`.
3. Fill three URLs from your cloud service:

- `submit_url`
- `status_url_template`
- `result_url_template`

4. (Optional) If your service needs auth, set env var:

```bash
set CLOUD_RF_API_KEY=your_token_here
```

5. Restart app.

## Check status

- `GET /api/cloud_status`
- `GET /api/rfdiffusion_status`

In UI, click `Generate Mini-binder` (Cloud), then click `Check Cloud Job` until completed.

## Expected cloud API contract

### Submit

`POST submit_url`

Request body:

```json
{
  "design_id": "PDCD1_20260402T120000Z",
  "receptor": "PDCD1",
  "cell": "CD8 T cell",
  "process": "exhaustion",
  "binder_length": 42,
  "hotspot_hint": "FG loop"
}
```

Response:

```json
{
  "job_id": "job_123",
  "status": "queued"
}
```

### Poll status

`GET status_url_template` with `{job_id}` replaced.

Response example:

```json
{
  "status": "running"
}
```

or

```json
{
  "status": "completed"
}
```

### Fetch result

`GET result_url_template` with `{job_id}` replaced.

Response should include design payload fields used by UI, especially:

- `design_id`
- `receptor`
- `engine`
- `binder_sequence`
- `binding`
- `files`
- `complex_pdb_text` (for immediate 3D rendering)

## Notes

- This app can still use mock output when cloud/local RFdiffusion is not available.
- Real RFdiffusion results require computational and experimental validation.
