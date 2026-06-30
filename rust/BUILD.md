# WASM 빌드 & docs/pkg_v5 교체 가이드

## 0. 사전 준비 (한 번만)
```bash
# Rust 툴체인
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source "$HOME/.cargo/env"
rustup target add wasm32-unknown-unknown
# wasm-pack
cargo install wasm-pack         # 또는: curl https://rustwasm.github.io/wasm-pack/installer/init.sh -sSf | sh
```

## 1. 빌드 차단 요인 — 이미 수정함 ✅
기존 코드는 `wasm-bindgen`의 제거된 기능(`serde-serialize` + `JsValue::from_serde`)을 써서
잠긴 `wasm-bindgen 0.2.100`으로는 **컴파일이 안 됐다.** 다음과 같이 고쳐둠:
- `Cargo.toml`: `wasm-bindgen = "0.2"` (serde-serialize 제거), `serde-wasm-bindgen = "0.6"` 추가
- `src/lib.rs`: 반환을 `result.serialize(&serde_wasm_bindgen::Serializer::json_compatible()).unwrap()` 로 변경
  - `json_compatible()` = Rust 맵을 JS **Object**로 직렬화 → 기존 JS(`Object.entries`, 속성 접근)가 그대로 동작

> 되돌리려면 백업: `work/Cargo.toml.bak`, `work/lib.rs.bak`

## 2. 빌드 (web 타깃 — 기존 glue가 import.meta/initSync를 쓰므로 반드시 web)
```bash
cd rust
wasm-pack build --release --target web --out-dir pkg_build
```
산출물(crate명 = CRISPRlungo_regular): `pkg_build/CRISPRlungo_regular.js`, `CRISPRlungo_regular_bg.wasm`, `*.d.ts`, `package.json`

## 3. docs/pkg_v5 교체
로더(`docs/lungo_worker.js`)는 `./pkg_v5/CRISPRlungo_regular.js` + `CRISPRlungo_regular_bg.wasm` 만 불러온다.
```bash
cp pkg_build/CRISPRlungo_regular.js          ../CRISPRlungo/docs/pkg_v5/
cp pkg_build/CRISPRlungo_regular_bg.wasm      ../CRISPRlungo/docs/pkg_v5/
cp pkg_build/CRISPRlungo_regular.d.ts         ../CRISPRlungo/docs/pkg_v5/      # 선택
cp pkg_build/CRISPRlungo_regular_bg.wasm.d.ts ../CRISPRlungo/docs/pkg_v5/      # 선택
```
- `CRISPRlungo_regular_new_*.{js,wasm}` 는 로더가 안 쓰는 구버전/중복 → 무시 또는 삭제.
- 한 번에 덮어쓰려면: `wasm-pack build --release --target web --out-dir ../CRISPRlungo/docs/pkg_v5`
  (단, 이러면 pkg_v5의 package.json 등도 덮어써짐 — 보통 무방)

## 4. 동작 확인
브라우저에서 `docs/index.html`을 (로컬 서버로) 열어 Example_Treated/Control 실행.
정적/로컬서버 예: `cd CRISPRlungo/docs && python3 -m http.server 8000` → http://localhost:8000

## 5. Python 파리티 재검증 (선택, 강력 추천)
재빌드한 wasm이 Python과 같은 결과를 내는지 자동 비교:
```bash
# (a) Python으로 같은 데이터 실행해 SAM + 기준 결과 생성
CRISPRlungo PD1.fasta --control Example_Control.fastq Example_Treated.fastq OUT ggcgccctggccagtcgtct
#    OUT/results/input_summary.txt 에서 CleavagePos 확인 (Example/PD1 = 2369, 단일 타깃)
# (b) wasm 구동 + 비교
node parity_check/run.mjs ../CRISPRlungo/docs/pkg_v5 OUT/align PD1.fasta 2369 none 5 /tmp/wasm_out.json
python3 parity_check/compare.py /tmp/wasm_out.json OUT/results
#    기대: identical_set=True, count mismatches: 0
```
직전(수정 전 wasm) 기준 결과: 유의집합 30=30 일치, read카운트 일치, 키메라 1건만 차이
(`2363:Del_19` vs `2368:Del_8`). R1 수정이 이 1건을 없애는지 위 비교로 확인.

## 트러블슈팅
- `error: failed to select a version ... serde-wasm-bindgen`: 인터넷(crates.io) 접속 필요. 사내망이면 프록시 설정.
- `feature 'serde-serialize' does not exist`: 1번 수정이 안 들어간 것 — Cargo.toml 확인.
- `getrandom` 관련 wasm 에러: 이미 `getrandom = { features=["js"] }` 설정돼 있어 정상.
- main.rs(네이티브 bin)는 `wasm-pack build --lib`에 포함 안 되니 무시해도 됨.
