# v4 最終投稿準備パッケージ

1. PPIをDeclarationsからMethods末尾へ移動。全過程の事実確認は著者入力欄として保持し、確認後に使える全文を`manuscript/submission_checklist.md`に用意。
2. Figure 3 Panel Bを既存B2結果で再描画。期待値曲線とmasking rootを元CSVに照合。青実線／赤破線を維持。補足図S4のB1/B2比較とS3の既存B1表示は変更していない。
3. `pytest.ini`に`pythonpath = .`を追加。consoleの`pytest -q`と`python -m pytest -q`の両方で16件合格。
4. 著者・所属、責任著者、資金、利益相反、貢献／guarantor、倫理判断、PPI、repository URL/DOI、AI model/version/dateと人手確認、cover letterを入力先付きチェックリストに整理。未確認事項は創作していない。
5. Cover letterの“subject to confirmation of the current article category”を削除。投稿画面の正式選択肢と字数・図数制限の最終確認はチェックリストに保持。
6. 新規simulationなし。全保存結果、全表、設定、図生成以外の解析source、protocol、Figure 3以外の全図をSHA-256で不変確認。

Song訂正は、提供されたv3最終監査に基づき解決済みと記録した（v4で原文を再確認したとは記載していない）。Word本文を更新し、全13ページとB2図を目視確認済み。現在の語数には著者向け確認欄を含むため、著者入力後に再計測する。

v1–v3のZIPは保持。過去auditは履歴。現版の検証は`finalisation_verification_v4.json`と`test_results_v4.txt`を参照。全体ZIPはcache・一時ファイル・描画中間物を含まない。投稿・公開depositは未実施。
