# AP Novo 입력 준비 (Rubisco bind_core / bind_c40 / bind_c45)

`inputs/` 의 세 구조를 AP Novo 입력용으로 정리한 것. Garibaldi에서 이 브랜치를 pull 해서 쓴다.

```bash
git fetch origin claude/alphaprotein-novo-analysis-3i26e2
git checkout claude/alphaprotein-novo-analysis-3i26e2
cd ap_novo_prep
python3 prep_ap_novo.py --selftest
python3 prep_ap_novo.py --inputs inputs/bind_core.pdb inputs/bind_c40.pdb inputs/bind_c45.pdb --out out
```

stdlib만 쓰므로 설치할 것이 없다. `out/` 은 이미 생성돼 있어서 다시 돌리지 않아도 된다.

## 하는 일
- 단백질 chain A/B는 원래 번호 그대로 유지
- KCX191 → LYS + FMT 리간드 (9RUB 방식). NZ–CX 결합(1.34 Å)은 입력에 표현되지 않는다
- CAP → chain C, MG → chain D, FMT → chain E (각 1 잔기)
- CONECT는 버린다. 리간드 토폴로지는 CCD 이름으로 해석된다고 가정
- `out/<이름>/` : `motif.pdb`, `motif_spec.json`, `report.txt`. `out/cross_check.txt` : 세 구조 비교

## 확인된 것 (이 컨테이너에서 실제로 실행)
- 입력 → 출력 좌표 보존: 단백질 원자 불일치 0, FMT 좌표 = 원래 카바메이트 CX/OQ1/OQ2
- 세 구조는 core ⊂ c40 ⊂ c45, 공유 원자 좌표 편차 0.000 Å
- 모든 잔기의 중원자 세트가 완전함
- Mg 6배위(2.03–2.17 Å): KCX OQ2, D193 OD1, E194 OE1, CAP O2/O3/O6
- 모티프 토큰 수: core 41, c40 45, c45 49 (모두 unindexed 상한 128 이내라고 앞선 메모에 있음)
- `--selftest`: T13/A10 레이아웃의 간격 계산이 앞선 분석의 motif_str을 재현

## 확인되지 않은 것
- **CAP 원자명 ↔ CCD**: RCSB/EBI/PDBj 모두 네트워크 정책으로 막혀 있어 대조 못 함. 이름은 O1P…O6P, O6, O7 형식
- **AP Novo manifest 문법**: README를 못 구해서 전체 manifest는 만들지 않았다. `motif_spec.json` 의
  `unverified_fragments_from_earlier_readme_reading` 은 앞선 세션이 README에서 읽었다는 필드(`motif_str`, `seq_length`,
  `unindexed_motif_residues`, `is_author_naming`)만 옮긴 조각이다. 키 이름과 값 형식은 실제 README로 확인할 것
- **indexed 레이아웃**: A10 번호표가 없어 만들지 못했다. `{"A164": 25, ...}` 형식의 native→설계 위치 JSON을
  `--numbering-map` 으로 주면 `motif_str` 이 생성된다 (모티프 잔기 전부가 들어 있어야 함)
- AP Novo 환경 설치, 가중치, 실제 실행: 하지 않았다

## 주의
- bind_c40/c45는 core에 없는 잔기(I164, M330, A392, E48' 등)를 포함한다. 잔기가 많을수록 어려워진다는
  앞선 분석 내용이 그대로 적용된다
- unindexed 에서는 출력의 모티프 좌표가 강제되지 않으므로 채점은 refold 구조 기준으로 해야 한다
