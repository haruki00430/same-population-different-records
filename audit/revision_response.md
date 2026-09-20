# 第2版：監査コメントへの対応

## 1. Introductionと比較対象

Introductionを4段落で再構成し、参考文献を12本に拡充しました。Abelらはrandom variationによる組織performanceの誤分類、Hoferは目的・条件・集団に依存するsignal/noise reliabilityとして位置付けています。本研究はこれらを否定せず、record-generation driftが経時的な改善判断を変えるという別の問いを扱います。

Sjodingらの先行recoding simulationも明示し、「記録変更のシミュレーション自体が新しい」という主張を避けました。SPC方法論とRhee、Songらの一次研究も追加・整理しました。文献ごとの確認範囲は `manuscript/citation_metadata.json` に保存しています。

## 2. Tipping pointを中心結果へ

Abstract、Results、Discussionの冒頭を、真の10%改善を期待値上消す8.889 ppの安全性ascertainment増加と7.273 ppのprocess numerator capture低下に変更しました。4.0%→4.5%は解析解の確認・説明例として残しています。有限標本のmasking確率と期待値のmasking rootは区別しています。

## 3. Parameter plausibility

Amaravadiらの出版社本文・Table 3とMedPACの一次資料を確認し、記録変化を感度分析で検討する動機としてDiscussionに追加しました。外部研究の周辺的な診断記録確率・重症度構成比は、本研究の真の疾患・事象を条件とする捕捉感度とは異なります。数値が近いことを8.889 ppの実証的検証とは解釈せず、較正したとも記述していません。政策時点も実施だけでなくannouncementを含むことを反映しました。

## 4. Risk adjustmentの位置付け

B1/B2の10%改善消失には14.658–18.447 ppを要し、安全性・プロセス指標より大きいことを明示しました。第3の指標構造への適用例として保持し、20%改善のrootがmain grid外であることも維持しています。

## 5. SPC comparatorと全応答

11.06%は24-month follow-up中に8 consecutive points ruleを少なくとも一度満たす確率です。1回の検定・1か月あたりの誤り確率とは区別し、Methods、Results、Discussion、Supplementを同期しました。

既存のSPC結果から補足図S5と表S9を作成し、80%のno-change comparatorと82%、85%、90%、95%の全endpointをabrupt/gradual・両ルールについて表示しました。重複する17個の8-point windowがあり、判定が独立ではないことも説明しました。新たなsimulation cellは追加していません。

## 6. タイトルと投稿先

タイトルを **Same Population, Different Records: Stress-Testing Healthcare Quality Measures Against Changes in Record Generation—a Simulation Study** に変更し、カバーレターも同期しました。BMJ Quality & Safety、Research and Reporting Methodologyを維持しています。記事種別は2026年Hofer論文のpublisher-indexed情報で確認しました。詳細なword/figure limit等は公式ページが取得できず、投稿直前の確認事項として残しています。

## 保持・検証

研究設計、設定、simulation source、結果summary、事前プロトコル、主要図1–4を保持し、保存済み結果の再表示と原稿生成のみを行いました。SHA-256照合、引用番号、補足表と元データの照合、既存16件のテストの結果は `revision_verification.json` と `test_results.txt` に保存しました。

元の解析auditは解析時点の記録として保持します。本改稿は解析後の編集上の改稿であり、事前登録を行ったとは主張していません。著者情報、Songの訂正情報、AHRQ発行年等の最終確認は残っています。

## 編集者への2文の回答

Previous studies establish that coding changes can alter specific measures and that random noise can misclassify organisational performance. This study supplies a common counterfactual stress-test across measure architectures that separates systematic record drift from clinical change and quantifies the recording shifts sufficient to erase or reverse a genuine improvement.
