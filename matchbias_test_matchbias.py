#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
matchbias 单元测试

运行：
    python -m unittest test_matchbias -v
    python -m pytest test_matchbias.py -v
"""
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from matchbias import anchored_pattern, unanchored_pattern, classify, DEFAULT_RULES


class TestPatterns(unittest.TestCase):
    """匹配模式的核心行为"""

    def test_anchored_rejects_substring(self):
        """词边界应拒绝跨基因子串：VEGFR 不该被 EGFR 命中"""
        pat = anchored_pattern("EGFR")
        self.assertFalse(pat.search("VEGFR"))
        self.assertFalse(pat.search("VEGFR1"))
        self.assertFalse(pat.search("VEGFR2"))
        self.assertFalse(pat.search("VEGFR3"))

    def test_unanchored_accepts_substring(self):
        """无词边界会误命中：VEGFR 被 EGFR 命中（这就是偏倚）"""
        pat = unanchored_pattern("EGFR")
        self.assertTrue(pat.search("VEGFR"))
        self.assertTrue(pat.search("VEGFR1"))

    def test_anchored_accepts_numeric_suffix(self):
        """词边界应保留数字后缀：EGFR2 仍是 EGFR"""
        pat = anchored_pattern("EGFR")
        self.assertTrue(pat.search("EGFR"))
        self.assertTrue(pat.search("EGFR2"))
        self.assertTrue(pat.search("EGFR3"))

    def test_anchored_rejects_multi_letter_suffix(self):
        """后缀是字母时视为不同基因（EGFRXY ≠ EGFR）"""
        pat = anchored_pattern("EGFR")
        self.assertFalse(pat.search("EGFRXY"))
        self.assertFalse(pat.search("EGFRB"))

    def test_case_insensitive_by_caller(self):
        """调用方应传入大写文本"""
        pat = anchored_pattern("EGFR")
        self.assertFalse(pat.search("egfr"))  # 小写不匹配（调用方负责 upper）


class TestClassify(unittest.TestCase):
    """归类逻辑"""

    def test_vegfrinhibitor_unanchored_goes_to_egfr(self):
        """无锚定：VEGFR 抑制剂被误归 EGFR/ERBB"""
        r = classify("PDGFR, KIT, VEGFR, RAF", DEFAULT_RULES, "Sorafenib", anchor=False)
        self.assertEqual(r, "EGFR/ERBB")

    def test_vegfrinhibitor_anchored_goes_to_other(self):
        """锚定：VEGFR 抑制剂正确归入 Other"""
        r = classify("PDGFR, KIT, VEGFR, RAF", DEFAULT_RULES, "Sorafenib", anchor=True)
        self.assertEqual(r, "Other")

    def test_true_egfr_stays_in_class(self):
        """真 EGFR 抑制剂两种方式都归 EGFR/ERBB"""
        for anchor in (False, True):
            r = classify("EGFR", DEFAULT_RULES, "Erlotinib", anchor=anchor)
            self.assertEqual(r, "EGFR/ERBB")

    def test_erk_is_meK_mapk(self):
        r = classify("ERK1, ERK2", DEFAULT_RULES, "Ulixertinib", anchor=True)
        self.assertEqual(r, "MEK/MAPK")

    def test_sl0101_special_case(self):
        """SL0101 通过特殊规则归入 AURK/mitotic"""
        r = classify("", DEFAULT_RULES, "SL0101", anchor=True)
        self.assertEqual(r, "AURK/mitotic")


class TestIntegration(unittest.TestCase):
    """端到端：用内置示例数据跑一遍"""

    def test_example_file_detects_bias(self):
        example = Path(__file__).parent / "example_gdsc1.csv"
        if not example.exists():
            self.skipTest("示例数据不存在")
        import subprocess
        r = subprocess.run(
            [sys.executable, str(Path(__file__).parent / "matchbias.py"),
             "--input", str(example), "--pathway-col", "Pathway"],
            capture_output=True, text=True
        )
        # 退出码 1 = 检出偏倚
        self.assertEqual(r.returncode, 1)
        self.assertIn("EGFR/ERBB", r.stdout)
        self.assertIn("Sunitinib", r.stdout)
        self.assertIn("Sorafenib", r.stdout)


if __name__ == "__main__":
    unittest.main(verbosity=2)
