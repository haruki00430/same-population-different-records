# 第3版の変更と確認

- 図1：旧版の矢印が箱の中から始まる・先端が箱に隠れる問題を修正。箱高、中心線、箱間隔を統一し、4本すべての矢印を箱の外側に配置。余白と注記位置を調整。PNG・SVG・PDFを更新し、PNGとWord配置を目視確認。
- Word本文：`manuscript/manuscript.docx`を追加。本文、抄録、参考文献、声明案、図説明、主要図1–4を収録。数式は編集可能なWord数式。全13ページをMicrosoft Wordで描画して確認。表・補足資料は従来どおり別ファイル。Word生成コードを同梱。
- 添付監査：Abstract/Discussionを8.9、7.3、14.7–18.4へ丸め、ResultsとTable 2に精密値を保持。B2を実務的主例、B1を機序のreferenceと明記。sentence-case titleに変更。図3を青実線／赤破線に変更。
- AI利用：実際の支援範囲をMethodsとContributorsに記載。著者による独立検証、PPI全過程、所属機関の倫理判断、著者情報、公開repositoryの実depositは確認欄として残した。未実施の公開・承認・人手検証を完了と記載していない。
- 新規simulationなし。設定、計算結果、解析コード（図生成以外）、全結果表、主要図2・4とprotocolはSHA-256で不変確認。既存16テストが合格。過去のv2 auditは履歴であり、現版の検証は`finalisation_verification_v3.json`を参照。
- 配布ZIPは一時ファイル・キャッシュ・レンダリング中間物を除いたlocal clean release。元のv1/v2 ZIPは保持。
