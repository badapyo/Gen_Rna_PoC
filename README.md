# Gen_Rna_PoC
AI 기반 RNA 생성 모델(Proof of Concept)

🚀 Overview

이 코드는 실제 RNA-seq 데이터(엑셀 파일)를 읽어와서,
1. Normal / Tumor 샘플을 짝(pair) 으로 묶고
2. 표현형을 잘 보존하도록 전처리 + 차원축소 를 한 뒤
3. Normal과 Tumor 각각에 대해 Bayesian Gaussian Mixture Model(BGM) 을 학습하고
4. 학습된 BGM으로부터 synthetic RNA 샘플(합성 Normal/Tumor) 을 생성해서
5. 다시 원래 gene expression 스케일(count 수준) 로 되돌린 후
6. 엑셀 파일(Synthetic_RNAseq_FullInverse.xlsx)로 저장
7. 마지막으로 시각화(Kernel PCA), MMD, 분류 성능(Logistic Regression) 으로
8. 합성 데이터가 실제 데이터와 얼마나 비슷한지 평가한다.


🧬 데이터 로딩 & Normal/Tumor 페어 구성

- 41467_2024_54434_MOESM11_ESM.xlsx 파일에서 Gene_ID 컬럼을 기준으로 RNA-seq matrix를 읽어온다.

- 샘플 컬럼명 규칙:
  - "<PAIR_ID>-N" : Normal
  - "<PAIR_ID>-T" : Tumor
    예) TCGA-AB-1234-N, TCGA-AB-1234-T

- 이걸 이용해서 각 열을
  - pair_id (같은 환자/케이스)
  - cond (N/T) 로 파싱해서

- Normal과 Tumor가 둘 다 존재하는 pair만 필터링 해서 사용한다.
- 최종적으로:
  - Xn : Normal matrix (pairs × genes)
  - Xt : Tumor matrix (pairs × genes)
    로 분리된다.


⚙️ 전처리 & 차원 축소

주요 하이퍼파라미터:
- TOP_VAR_GENES = 8000 : 변동이 큰 상위 유전자만 사용
- USE_PCA = True
- PCA_DIM = 128
- SEED = 42

전처리 흐름:

1. 숫자형 변환(to_num)
  - 문자열 섞여 있는 값을 float로 강제 변환, 실패 시 NaN 부여.
2. log1p 변환
  - log1p(count) 로 스케일을 안정화.
3. StandardScaler
  - 각 gene에 대해 평균 0, 분산 1로 표준화.
4. 상위 변동 유전자 선택
  - Normal/Tumor 전체에서 variance가 큰 상위 TOP_VAR_GENES만 남김.
5. PCA (선택적)
  -  PCA_DIM 차원까지 축소해서 Xn_m, Xt_m (model/latent space) 생성.
  -  PCA를 사용하지 않으면, 표준화된 gene space 그대로 사용.


📦 Bayesian Gaussian Mixture Model(BGM) 학습
- fit_bgm(X) 함수에서 BayesianGaussianMixture 모델을 학습.
- Normal/Tumor 각각에 대해 별도로 학습:
  - bgm_N = fit_bgm(Xn_m)
  - bgm_T = fit_bgm(Xt_m)
- BGM 특징:
  - Dirichlet Process prior 때문에 필요 없는 component weight는 0 근처로 수렴.
  - 최대 component 개수보다 실제로는 적은 수의 active component 사용.
  - 일반 GMM보다 과적합이 덜하고, 작은 데이터에도 안정적.



🎲 합성 샘플 생성 (Synthetic RNA)
- 클래스별 합성 샘플 수: N_SYN = 500
- sample_bgm(model, n):
  - BGM으로부터 n개의 샘플을 샘플링.
- Normal/Tumor 합성 latent 샘플:
  - synN = sample_bgm(bgm_N, N_SYN)
  - synT = sample_bgm(bgm_T, N_SYN)

✅ 분포 보정 (moment / covariance matching)
- psd_sqrt(mat):
  - 대칭 양의 준정정부호(PSD) 행렬의 matrix square root를 SVD로 계산.
- (코드 안의 match_moments 류 함수):
  - 합성 샘플의 평균·공분산 구조를 실제 데이터와 맞추도록 affine transform 수행.
  - Real Normal/Tumor의 1·2차 통계 구조를 더 잘 따라가도록 조정.



🔁 latent space → gene expression 스케일로 역변환
- inverse_to_gene(Xm) / inverse_to_original(X_model) 함수가 핵심.

처리 순서:
1. PCA 역변환
  - latent/model space → PCA 이전의 batch-scaled gene space
2. StandardScaler 역변환
  - 표준화 해제 → log1p gene space로 복원
3. np.expm1
  - log1p(count) → count 스케일로 복원
4. NaN / 음수 처리
  - NaN 은 0으로 치환
  - 음수는 0으로 클리핑 (RNA count는 음수가 될 수 없음)

결과:
- synN_gene, synT_gene : (samples × genes)의 합성 gene expression matrix.



📊 결과 저장 (Synthetic RNA Excel)
- 합성 Normal/Tumor 데이터를 DataFrame으로 변환:
  - index: Gene_ID
  - columns: SynN_0001, SynN_0002, … / SynT_0001, …
- Gene_ID 를 첫 컬럼으로 옮겨서 원본 포맷과 최대한 비슷하게 정리.



👀 시각화: Kernel PCA 2D embedding
- KernelPCA(RBF kernel)를 이용해서
  - Real Normal (RN),
  - Real Tumor (RT),
  - Synthetic Normal (SN),
  - Synthetic Tumor (ST)
    를 한 번에 2D로 투영.
- auto_gamma(X):
  - RBF kernel의 gamma 를 데이터의 pairwise 거리 중앙값(median heuristic)으로 자동 추정.
- 2D scatter plot에서
  - Real vs Synthetic 분포가 얼마나 겹치는지 시각적으로 확인.



📐 통계적 거리: MMD (Maximum Mean Discrepancy)
- mmd_rbf(X, Y, gamma=None):
  - RBF kernel 기반 MMD² 계산.
  - rbf_kernel 로 Gram matrix를 만들고,
  - E[XX] + E[YY] - 2E[XY] 형태로 MMD² 구함.
- Real vs Synthetic에 대해 각각 계산:
  - mmdN, gN = mmd_rbf(RN, SN)
  - mmdT, gT = mmd_rbf(RT, ST)


🧪 분류 기반 평가: Real vs Synthetic 구분 성능
- eval_cross(Xtr, ytr, Xte, yte, name):
  - StandardScaler + LogisticRegression 파이프라인으로 이진 분류.
  - roc_auc_score, f1_score, accuracy_score 출력.
- 대표적인 실험:
  - Real로 학습 → Synthetic에 테스트
  - Synthetic으로 학습 → Real에 테스트
- 만약 분류기가 Real vs Synthetic을 잘 구분 못하면(낮은 AUC/F1)
  → Synthetic 데이터가 Real 데이터와 통계적으로 비슷하다는 신호로 해석 가능.


